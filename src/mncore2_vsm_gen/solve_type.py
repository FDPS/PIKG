import copy
from lark import Lark, Transformer

def solve_type(instlist,decllist, verbose = False):
    tdic = {}
    matchtypes = MatchTypes(TYPERULES) # Defined Below
    for decl in decllist:
        TYPE = decl[1][0][1]
        name = decl[1][1][1]
        TYPE = TYPE + "3" if TYPE[-3:] == "vec" else TYPE 
        tdic["M#"+name] = [TYPE]
    for op, ili, oli in instlist:
        for name in ili+oli:
            if name not in tdic:
                tdic[name] = None
    print ( tdic ) if verbose else None
    for loop_count in range(65536):
        tdic_old = copy.deepcopy(tdic)
        for op, ili, oli in instlist:
            #print ("--",op,ili,"->",oli)
            matchtypes.solve( [op, ili, [oli[0]] if len(oli)>0 else []], tdic ) # assuming that len(oli)>1 happens only when writes to mem and mask
        if tdic == tdic_old:
            break
        print ( tdic ) if verbose else None
    if all([ (v is not None) and (len(v)==1) for v in tdic.values()]):
        print ("type successfully solved") if verbose else None
    else:
        assert 1==0
        print ("DAMEDESITA")
    new_instlist = []
    for op, ili, oli in instlist:
        new_instlist.append( (op,[i+":"+tdic[i][0] for i in ili] ,[o+":"+tdic[o][0] for o in oli]) )
    return new_instlist

class MatchTypes:
    def __init__(me,TYPERULES):
        checker_parser = Lark(TYPECHECKER_BNF,start='program')
        tree = checker_parser.parse(TYPERULES)
        tree = enTree().transform(tree)
        me.rules = expandRules(tree)
    def solve(me,inst, tdic):
        op, ili, oli = inst
        matched_rules = []
        for rule in me.rules:
            r_op, r_ili, r_oli = rule
            if op != r_op:
                continue
            if len(r_ili) != len(ili):
                continue
            if len(r_oli) != len(oli):
                continue
            ilichk = [ (tdic[i] is None) or (ri in tdic[i]) for i,ri in zip(ili,r_ili) ]
            olichk = [ (tdic[o] is None) or (ro in tdic[o]) for o,ro in zip(oli,r_oli) ]
            if all( ilichk + olichk ):
                matched_rules.append(rule)
        if len(matched_rules) == 0:
            print (tdic)
            print ("Err", inst )
            assert 0=="Type inference failed. There are missing rules or inconsistencies."
        else:
            extdic={ it:set() for it in ili+oli}
            for rule in matched_rules:
                r_op, r_ili, r_oli = rule
                for it,tp in zip(ili+oli,r_ili+r_oli):
                    extdic[it].add(tp)
            for k,v in extdic.items():
                tdic[k] = list(v)
                
class enTree(Transformer):
    def __default__(me,data,children,meta):
        return(data,children)
    def __default_token__(me,token):
        return (token.value)
    
from itertools import product
def expandRules(tree):
    tempdec = {}
    rules = []
    instlist = tree[1]
    for head, container in instlist:
        if head == "tempdec":
            tempdec[container[0]] = container[1:]
            #print (tempdec)
        if head == "type_statement":
            (_,tlist), fname, (_,ili), (_,oli) = container
            tlist = [ b for a,b in tlist ]
            xrulelist = [ (fname, ili, oli) ]
            for template, from_name in tlist:
                nrulelist = []
                to_name_li = tempdec[template]
                for to_name, (fname, ili, oli) in product(to_name_li,xrulelist) :
                    ili = [ to_name if it==from_name else it for it in ili ]
                    oli = [ to_name if it==from_name else it  for it in oli ]
                    nrulelist.append( (fname,ili,oli) )
                xrulelist = nrulelist
            rules.extend( xrulelist )
    return rules

TYPECHECKER_BNF = """
    program : (statement? _NEWLINE)* statement? -> program
    
    ?statement: type_statement
              | tempdec_statement
    
    type_statement : template_list FNAME ":" typelist "->" typelist -> type_statement
    template_list : template* -> template_list
    template : "<" TYPE ":" TYPE ">" -> template
    typelist: (TYPE? ",")* TYPE? -> typelist
    
    tempdec_statement: "template" TYPE "=" TYPE ("," TYPE)* -> tempdec
    
    TYPE : CNAME
    FNAME : CNAME
    
    
    _NEWLINE : NEWLINE
    _WS_INLINE : WS_INLINE

    %import common.CNAME
    %import common.NEWLINE
    %import common.WS_INLINE
    %ignore WS_INLINE
"""

TYPERULES = """
template IFFv = S16, S32, S64, F16, F32, F64, F16vec2, F32vec2, F64vec2, F16vec3, F32vec3, F64vec3, F16vec4, F32vec4, F64vec4
template FFv = F16, F32, F64, F16vec2, F32vec2, F64vec2, F16vec3, F32vec3, F64vec3, F16vec4, F32vec4, F64vec4
template IF = S16, S32, S64, F16, F32, F64
template F = F16, F32, F64
template F16v = F16vec2, F16vec3, F16vec4
template F32v = F32vec2, F32vec3, F32vec4
template F64v = F64vec2, F64vec3, F64vec4
template Fv = F16vec2, F16vec3, F16vec4, F32vec2, F32vec3, F32vec4, F64vec2, F64vec3, F64vec4

<IFFv:type> add: type, type -> type
<IFFv:type> sub: type, type -> type
<IFFv:type> set: type -> type

<IF:type> mul: type, type -> type
<F16v:type> mul: type, type -> F16
<F32v:type> mul: type, type -> F32
<F64v:type> mul: type, type -> F64
<F16v:type> mul: type, F16 -> type
<F32v:type> mul: type, F32 -> type
<F64v:type> mul: type, F64 -> type
<F16v:type> mul: F16, type -> type
<F32v:type> mul: F32, type -> type
<F64v:type> mul: F64, type -> type

<F:type> rsqrt: type->type

<F:type>    to_f32:  type->F32
<Fv:type> to_f32vec: type->F32vec3
<F:type>    to_f64:  type->F64
<Fv:type> to_f64vec: type->F64vec3

<IF:type> imm: type -> type
<IF:type> cmp: type, type -> type

<IF:type> and: type, type -> type
<IF:type> xor: type, type -> type
<IF:type> dec: type -> type

lnot: BOOL -> BOOL
<IF:type> lor: type, type -> BOOL
<IF:type> lor: BOOL, type -> BOOL
<IF:type> lor: type, BOOL -> BOOL
<IF:type> land: type, type -> BOOL
<IF:type> land: BOOl, type -> BOOL
<IF:type> land: type, BOOL -> BOOL


mask: BOOL ->
mend: ->
"""
