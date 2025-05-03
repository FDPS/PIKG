import sys
import json

from add_build_in_functions import *
from code_parser            import *
from add_decllist_from_json import *
from remove_redundant_set   import *
from solve_type             import *
from remove_type_conversion import *
from generate_core_instlist import *
from code_annealer          import *
from generate_head_instlist import *

filename = sys.argv[1]

with open(filename) as f:
    main_code = f.read() 
code = add_build_in_functions(main_code)
tree_trace_result = code_parser(code)
tree_trace_result.decllist = add_decllist_from_json(tree_trace_result.decllist)
nlist = order_graph().run(tree_trace_result.instlist, tree_trace_result.decllist)
ntlist = solve_type(nlist,tree_trace_result.decllist)
ntlist = remove_type_conversion(ntlist)
core_instlist_gen, core_decllist, Rnum, Ilen, imm_instlist_gen = generate_core_instlist(ntlist, tree_trace_result.decllist)
core_instlist_al, found_state = code_annealer(core_instlist_gen, Rnum, Ilen, nloops=5)
head_obj = generate_head_instlist(core_instlist_al, found_state, imm_instlist_gen,core_decllist)

with open("codehead.json","w") as f:
    jsonstr = json.dumps(head_obj)
    f.write( jsonstr )
