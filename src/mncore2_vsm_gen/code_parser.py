from lark import Lark, Transformer

def code_parser(code):
    parsed_code = Lark(PIKG_BNF,start='program').parse(code)
    tree = makeTree().transform(parsed_code)
    tree_trace_result = treeTraceResult()
    tree_trace_result.sim(tree)
    return tree_trace_result

####
class treeTraceResult:
    def __init__(me):
        me.instlist = []
        me.decllist = []
        me.functions = {}
        me.fstack = ["M"]
        me.pstack = []
        me.funcid = 0
        me.pushid = 0
        me.maskid = 0
    def function_enter(me):
        me.funcid += 1
        fname = "F"+str(me.funcid)
        me.fstack.append(fname)
        return fname
    def function_back(me):
        me.fstack.pop()
    def mask_enter(me):
        me.maskid += 1
        return me.get_mask()
    def mask_back(me):
        me.maskid -= 1
    def get_mask(me):
        if me.maskid > 0:
            return "O#"+str(me.maskid)
        else:
            return None
    def var(me,cname):
        return me.fstack[-1] + "#" + cname
    def push(me):
        me.pushid +=1
        me.pstack.append(me.pushid)
        return "L#"+str(me.pstack[-1])
    def push_wo_inc(me):
        me.pstack.append(me.pushid)
        return "L#"+str(me.pstack[-1])
    def pop(me):
        return "L#"+str(me.pstack.pop())
    def run(me,op,ili,oli):
        ili = ili if isinstance(ili,list) else [ili]
        oli = oli if isinstance(oli,list) else [oli]
        ##print ("#",op,ili,"->",oli)
        me.instlist.append( [op,ili,oli] )
    def sim(me,tree):
        if tree is None:
            return
        
        name, child = tree
        
        if isinstance(child, list):
            if name == "function_block":
                if len(child) == 4:
                    funcname, varlist, progblock, ret_statement = child
                    me.functions[funcname[1]] = (varlist, progblock,ret_statement)
                elif len(child) == 3:
                    funcname, varlist, ret_statement = child
                    me.functions[funcname[1]] = (varlist, None,ret_statement)
                else:
                    assert 1=="function block statement illegal"
                return
            decl_names = ["decl_local", "decl_epi", "decl_epj", "decl_force", "decl_local_epi", "decl_local_epj","decl_local_force"]
            if name in decl_names:
                me.decllist.append(tree)
                return
            if name == "if_block":
                comp = child[0]
                blockTrue = None
                blockFalse = None
                for ch in child:
                    if ch[0] == "progblock":
                        blockTrue = ch
                    if ch[0] == "else_block":
                        blockFalse = ch
                # evaluate blockFalse first for easy masking
                if blockFalse is not None:
                    me.sim(blockFalse)
                me.mask_enter()
                me.sim(comp)
                flagvar = me.pop()
                if blockTrue is not None:
                    me.run("mask",me.get_mask(),[])
                    me.sim(blockTrue)
                    me.run("mend",[],[])
                me.mask_back()
                return
            
            # otherwise (evaluate child first)
            if name == "logical_and":
                me.sim(child[0])
                me.run('mask',[me.get_mask()],[])
                me.push_wo_inc()
                me.sim(child[1])
                me.run('mend',[],[])
                return
            else:
                for ch in child:
                    me.sim(ch)
            
            if name == "progblock" or name == "else_block" or name == "program": # this head is not meaningful (but child is)
                return
            if name == "return":
                return
            if name == "comp_rightnotlarge":  #"A >= B" === "cmp A B"
                reads = [me.pop(),me.pop()][::-1]
                writes = [me.push()]
                writes.append(me.get_mask())
                me.run( "cmp", reads, writes)
                return
            if name == "comp_leftnotlarge":  #"A <= B" === "cmp B A"
                reads = [me.pop(),me.pop()][::-1][::-1]
                writes = [me.push()]
                writes.append(me.get_mask())
                me.run( "cmp", reads, writes ) # reversed order
                return
            if name == "comp_leftlarge": #"A > B" === "static_cast<int>(A - B) - 1"
                me.run( "cmp", [me.pop(),me.pop()][::-1], me.push())
                reads = [me.pop()]
                writes = [me.push(), me.get_mask()]
                me.run( "dec", reads, writes)
                return
            if name == "comp_rightlarge": # A < B ==  static_cast<int>(B-A) - 1
                me.run( "cmp", [me.pop(),me.pop()][::-1][::-1], me.push())
                reads = [me.pop()]
                writes = [me.push(), me.get_mask()]
                me.run( "dec", reads, writes)
                return
            if name == "comp_equal":
                reads = [me.pop(),me.pop()]
                writes = [me.push()]
                writes.append(me.get_mask())
                me.run( "xor", reads, writes)
                return
            if name == "logical_or":
                print("logical_or(||) is not supported now")
                raise AssertionError()
                me.run( "lor", [me.pop(),me.pop()][::-1], me.push())
                return
            #if name == "logical_and":
            #    me.run( "land", [me.pop(),me.pop()][::-1], me.push())
            #    return
            if name == "assign":
                target_var = child[0][1]
                me.run( "set", me.pop(), me.var(target_var) )
                return
            if name == "assign_add":
                target_var = child[0][1]
                me.run( "add", [me.var(target_var),me.pop()], me.var(target_var) )
                return
            if name == "assign_sub":
                target_var = child[0][1]
                me.run( "sub", [me.var(target_var),me.pop()], me.var(target_var) )
                return
            if name == "add_op":
                me.run( "add", [me.pop(),me.pop()][::-1], me.push() )
                return
            if name == "sub_op":
                me.run( "sub", [me.pop(),me.pop()][::-1], me.push() )
                return
            if name == "mul_op":
                me.run( "mul", [me.pop(),me.pop()][::-1], me.push() )
                return            
            if name == "function_call":
                callname = child[0][1]
                if callname == "rsqrt_approx":
                    me.run( "rsqrt", me.pop(), me.push() ) ## approx code for GPFN2, thus newton process is requred.
                    return
                cast_functions=["to_f32","to_f32vec","to_f64","to_f64vec"]
                if callname in cast_functions:
                    me.run( callname, me.pop(), me.push() )
                    return
                else:
                    varlist, progblock, ret_statement = me.functions[callname]
                    fid = me.function_enter()
                    ##print ("--call", callname, fid)
                    for var in varlist[1][::-1]: # push reverse order from pop
                        target_var = var[1]
                        me.run( "set", me.pop(), me.var(target_var) )
                    me.sim( progblock )
                    me.sim( ret_statement )
                    me.function_back()
                    ##print ("--back", callname)
                    return
                print ("??call defined but not implemented?? "+callname)
            print ("??undefined operation??",name)
        else:
            if name == "WRITE_VAR":
                return
            if name == "READ_VAR":
                me.run( "set", me.var(child), me.push() )
                return
            if name == "VALUE_INT":
                me.run( "imm", f"VI#{child}", me.push() )
                return
            if name == "VALUE_FLOAT":
                me.run( "imm", f"VF#{child}", me.push() )
                return
            if name == "FUNCTION_NAME":
                return
            print ("??undefined leaf??",name,child)

#### makeTree class
class makeTree(Transformer):
    def __default__(me,data,children,meta):
        return(data,children)
    def __default_token__(me,token):
        return (token.type, token.value)

#### BNFs
PIKG_BNF = r"""
    ?program: (rootblock? _NEWLINE)* rootblock? -> program
    
    ?rootblock: declare_statement
              | function_block
              | progblock
              
              
    declare_statement: "EPI"   _WS_INLINE _type_and_var ":" MEMBER_VAR -> decl_epi
                     | "EPJ"   _WS_INLINE _type_and_var ":" MEMBER_VAR -> decl_epj
                     | "FORCE" _WS_INLINE _type_and_var ":" MEMBER_VAR -> decl_force
                     | "EPI"   _WS_INLINE "local" _WS_INLINE _type_and_var -> decl_local_epi
                     | "EPJ"   _WS_INLINE "local" _WS_INLINE _type_and_var -> decl_local_epj
                     | "FORCE" _WS_INLINE "local" _WS_INLINE _type_and_var -> decl_local_force
                     | _type_and_var                        -> decl_local
    _type_and_var: TYPEIS _WS_INLINE DECL_VAR
    
    
    function_block: "function" FUNCTION_NAME "("  function_varlist  ")" (_NEWLINE progblock)? _NEWLINE return_statement _NEWLINE* "end" -> function_block
    function_varlist: (FUNCTION_VAR ( "," FUNCTION_VAR )*)? -> function_varlist
    
    ?return_statement: "return" _WS_INLINE expression -> return
    
    
    ?progblock: (inner_statement? _NEWLINE)* inner_statement? -> progblock
    ?inner_statement: if_statement
                    | assign_statement

    ?if_statement: "if" if_block "endif"
    if_block: expression_bool (_NEWLINE progblock)? (_NEWLINE elseif_statement)? -> if_block
    ?elseif_statement: "else" (_NEWLINE progblock)? -> else_block
                    | "elseif" if_block -> else_block

    assign_statement: WRITE_VAR "=" expression -> assign
                    | WRITE_VAR "+=" expression -> assign_add
                    | WRITE_VAR "-=" expression -> assign_sub
                    | WRITE_VAR "*=" expression -> assign_mul
    
    ?expression_bool: expression_bool_or
    ?expression_bool_or: expression_bool_or "||" expression_bool_and -> logical_or
                       | expression_bool_and
    ?expression_bool_and: expression_bool_and "&&" expression_bool_comp -> logical_and
                        | expression_bool_comp
    ?expression_bool_comp: expression "==" expression -> comp_equal
                         | expression ">" expression -> comp_leftlarge
                         | expression "<" expression -> comp_rightlarge
                         | expression ">=" expression -> comp_rightnotlarge
                         | expression "<=" expression -> comp_leftnotlarge
                         | "(" expression_bool ")"

    
    ?expression: expression "+" term -> add_op
               | expression "-" term -> sub_op
               | term
    ?term: term "*" factor -> mul_op
         | factor
    ?factor: READ_VAR
           | VALUE_INT
           | VALUE_FLOAT
           | "(" expression ")"
           | FUNCTION_NAME "(" (expression ( "," expression )*)? ")" -> function_call
    
    
    TYPEIS: /(F(16|32|64)(vec(2|3|4)?)?)|((S|U)(16|32|64))/ 
    
    _NEWLINE: NEWLINE
    _WS_INLINE: WS_INLINE
    
    WRITE_VAR: CNAME
             | "FORCE" "." CNAME
    READ_VAR:  CNAME
             | "EPI" "." CNAME
             | "EPJ" "." CNAME
    MEMBER_VAR: CNAME
    DECL_VAR: CNAME
    FUNCTION_VAR: READ_VAR
    FUNCTION_NAME: CNAME
    
    VALUE_INT: SIGNED_INT
    VALUE_FLOAT: SIGNED_FLOAT
               | SIGNED_FLOAT "f"

    %import common.SIGNED_INT
    %import common.SIGNED_FLOAT
    %import common.CNAME
    %import common.NEWLINE
    %import common.WS_INLINE
    %ignore WS_INLINE

    """
