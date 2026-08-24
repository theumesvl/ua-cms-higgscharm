import json
import math
from pathlib import Path

import numpy as np


class TMVAGradBDT:
    """
    ROOT-free evaluator for TMVA BDT::BDTG models using BoostType=Grad.

    The implementation follows TMVA's GradBoost evaluation:

        sum = sum(tree_response)

        MVA = 2 / (1 + exp(-2 * sum)) - 1

    For GradBoost, the boostWeight stored in the XML is not used
    during evaluation. TMVA sums the response of the reached leaf
    from every tree. See line 1420-1429 in https://root.cern/doc/v632/MethodBDT_8cxx_source.html
    """

    def __init__(self, json_file):

        json_file = Path(json_file)

        with open(json_file, "r") as f:
            model = json.load(f)

        self.model = model

        self.ntrees = model["ntrees"]

        self.variables = model["variables"]

        self.variable_names = [
            variable["label"]
            for variable in self.variables
        ]

        self.trees = [
            tree["root"]
            for tree in model["trees"]
        ]

        # Sanity checks
        if len(self.trees) != self.ntrees:
            raise RuntimeError(
                f"Expected {self.ntrees} trees, "
                f"but found {len(self.trees)}"
            )

        boost_type = model["options"].get("BoostType")

        if boost_type != "Grad":
            raise RuntimeError(
                f"This evaluator expects BoostType=Grad, "
                f"but the XML contains BoostType={boost_type}"
            )

    # ------------------------------------------------------------
    # Single-tree evaluation
    # ------------------------------------------------------------

    @staticmethod
    def _evaluate_tree(node, values):

        while node["ivar"] >= 0:

            value = values[node["ivar"]]

            cut = node["cut"]

            # TMVA cType=1 corresponds to:
            #
            #     variable > cut
            #
            # going to the right child.
            #
            # Otherwise go left.
            if value > cut:
                node = node["right"]
            else:
                node = node["left"]

        return node["response"]

    # ------------------------------------------------------------
    # Single event evaluation
    # ------------------------------------------------------------

    def evaluate(self, values):

        """
        Evaluate one event.

        Parameters
        ----------
        values:
            Sequence of the 13 input variables in XML order.

        Returns
        -------
        float
            TMVA BDTG output in [-1, 1].
        """

        values = np.asarray(values, dtype=np.float64)

        if len(values) != len(self.variable_names):
            raise ValueError(
                f"Expected {len(self.variable_names)} variables, "
                f"got {len(values)}"
            )

        total = 0.0

        for tree in self.trees:
            total += self._evaluate_tree(tree, values)

        # TMVA:
        #
        # 2/(1+exp(-2*sum))-1
        #
        # Use tanh(sum), which is mathematically equivalent due to tanh(-x) = -tanh(x) and tanh(x) = (exp(2x)-1)/(exp(2x)+1)
        return math.tanh(total)

    # ------------------------------------------------------------
    # Evaluate many events
    # ------------------------------------------------------------

    def evaluate_array(self, values):

        """
        Evaluate an array of events.

        Parameters
        ----------
        values:
            numpy array with shape:

                (Nevents, 13)

        Returns
        -------
        numpy.ndarray
            Shape (Nevents,)
        """

        values = np.asarray(values, dtype=np.float64)

        if values.ndim != 2:
            raise ValueError(
                "values must be a 2D array"
            )

        if values.shape[1] != len(self.variable_names):
            raise ValueError(
                f"Expected {len(self.variable_names)} variables, "
                f"got {values.shape[1]}"
            )

        output = np.empty(
            values.shape[0],
            dtype=np.float32,
        )

        for i in range(values.shape[0]):
            output[i] = self.evaluate(values[i])

        return output
