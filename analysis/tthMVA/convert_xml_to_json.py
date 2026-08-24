import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def parse_xml(xml_file):
    xml_file = Path(xml_file)

    tree = ET.parse(xml_file)
    root = tree.getroot()

    # ------------------------------------------------------------
    # Variables
    # ------------------------------------------------------------

    variables_element = root.find("Variables")
    print("Variables element:")
    print(variables_element)

    variables = []

    for variable in variables_element.findall("Variable"):
        variables.append({
            "index": int(variable.attrib["VarIndex"]),
            "expression": variable.attrib["Expression"],
            "label": variable.attrib["Label"],
        })

    variables.sort(key=lambda x: x["index"])
    print("Variables:")
    print(variables)

    # ------------------------------------------------------------
    # Options
    # ------------------------------------------------------------

    options = {}

    options_element = root.find("Options")

    for option in options_element.findall("Option"):
        name = option.attrib["name"]
        value = option.text

        if value is None:
            continue

        value = value.strip()

        # Convert numerical options
        try:
            if "." in value or "e" in value.lower():
                value = float(value)
            else:
                value = int(value)
        except ValueError:
            pass

        options[name] = value

    # ------------------------------------------------------------
    # Trees
    # ------------------------------------------------------------

    weights_element = root.find("Weights")

    ntrees = int(weights_element.attrib["NTrees"])

    trees = []

    for binary_tree in weights_element.findall("BinaryTree"):

        tree_id = int(binary_tree.attrib["itree"])
        boost_weight = float(binary_tree.attrib["boostWeight"])

        root_node = binary_tree.find("Node")

        tree_data = {
            "itree": tree_id,
            "boostWeight": boost_weight,
            "root": parse_node(root_node),
        }

        trees.append(tree_data)

    trees.sort(key=lambda x: x["itree"])

    # ------------------------------------------------------------
    # Build output
    # ------------------------------------------------------------

    result = {
        "format": "TMVA_BDT_JSON",
        "method": root.attrib.get("Method", ""),
        "ntrees": ntrees,
        "variables": variables,
        "options": options,
        "trees": trees,
    }

    return result


def parse_node(node):

    result = {
        "ivar": int(node.attrib["IVar"]),
        "cut": float(node.attrib["Cut"]),
        "ctype": int(node.attrib["cType"]),
        "response": float(node.attrib["res"]),
        "ntype": int(node.attrib["nType"]),
    }

    children = list(node)

    if children:
        result["left"] = parse_node(children[0])
        result["right"] = parse_node(children[1])

    return result


def main():
    
    if len(sys.argv) != 3:
        print(
            "Usage:\n"
            "  python convert_xml_to_json.py input.xml output.json\n"
            "Or:\n"
            "  python3 convert_xml_to_json.py input.xml output.json"
        )
        sys.exit(1)
    cwd = Path.cwd()
    if cwd.name == "tthMVA":
        base_path = cwd.parent.parent
    elif cwd.name == "higgscharm" or cwd.name == "ua-cms-higgscharm":
        base_path = cwd.name
    else:
        raise ValueError(
            f"You are running the command 'convert_xml_to_json' from the wrong directory. Run it either from the main higgscharm dir or from analysis/tthMVA."
        )
    xml_file = (
        base_path
        / "analysis/data/tthMVA_2022-2023_retrained"
        / sys.argv[1]
    )

    print("xml_file path:")
    print(xml_file)

    json_file = (
        base_path
        / "analysis/data/tthMVA_2022-2023_retrained"
        / sys.argv[2]
    )
    print("json_file path:")
    print(json_file)
    result = parse_xml(xml_file)

    json_file.parent.mkdir(parents=True, exist_ok=True)

    with open(json_file, "w") as f:
        json.dump(result, f, separators=(",", ":"))

    print(f"Wrote {json_file}")
    print(f"Trees: {result['ntrees']}")
    print(f"Variables: {len(result['variables'])}")


if __name__ == "__main__":
    main()
