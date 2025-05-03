import json
def add_decllist_from_json(decllist):
    with open("decllist.json","r") as f:
        jsonstr = f.read()
        load_obj = json.loads(jsonstr)

        decllist_add = load_obj["decllist"]
        return decllist + decllist_add
