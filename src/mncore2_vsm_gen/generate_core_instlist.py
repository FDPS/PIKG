def generate_core_instlist(instlist,decllist):
    core_instlist = remove_local_var_calculation(instlist,decllist)
    core_decllist = remove_unused_var(core_instlist,decllist)
    core_instlist, imm_instlist = imm_divide(core_instlist)
    core_instlist_fma = fma_solver(core_instlist)
    core_instlist_uv = unvectorizeInst(core_instlist_fma)
    core_instlist_gen, Rnum, Ilen, imm_instlist_gen = generate_instlist_with_undetermined_raw_registers(
        core_instlist_uv, core_decllist, imm_instlist)
    core_instlist_gen = split_f64_for_f32(core_instlist_gen)
    core_instlist_gen, Rnum, Ilen, imm_instlist_gen = add_rotate_for_simd(core_instlist_gen,core_decllist,imm_instlist_gen, Rnum, Ilen)
    return core_instlist_gen, core_decllist, Rnum, Ilen, imm_instlist_gen
        
############################
def remove_unused_var(instlist,decllist):
    ret = []
    for decl in decllist:
        reg = decl[1][1][1]
        isUsed = False
        for op, reads, writes in instlist:
            for r in reads+writes:
                if r.split("#")[1].split(":")[0] == reg:
                    isUsed = True
                    break
        if isUsed:
            ret.append(decl) 
    return ret
def remove_local_var_calculation(instlist,decllist):
    #print(instlist)
    local_ep = []
    local_force = []
    force = []
    for decl in decllist:
       name = decl[1][1][1]
       if decl[0] in ['decl_local_epi','decl_local_epj']:
           local_ep.append(name) if name not in local_ep else None
       if decl[0] in ['decl_local_force']:
           local_force.append(name) if name not in local_force else None
       if decl[0] in ['decl_force']:
           force.append(name) if name not in force else None
    #print("local_ep",local_ep)
    #print("local_force",local_force)
    ret = []
    for inst in instlist:
        def get_name(item):
            return item.split("#")[1].split(":")[0]
        local_ep_computation = any([(get_name(w) in local_ep) for w in inst[2]])
        local_force_computation = False
        for r in inst[1]:
            for w in inst[2]:
                if (get_name(r) in local_force) and (get_name(w) in force):
                    local_force_computation = True
        if not (local_ep_computation or local_force_computation):
            ret.append(inst) 
        #else:
            #print("removed",inst)
    #print("ret",ret)
    return ret

def imm_divide(instlist):
    imm_instlist = []
    res_instlist = []
    for lno, (op, ili, oli) in enumerate(instlist):
        if "imm" in op:
            canDivide = True
            for i in range(lno+1,len(instlist)):
                _, _, oli_ = instlist[i]
                if oli_ == oli:
                    canDivide = False
            if canDivide:
                imm_instlist.append( [op, ili, oli] )
                continue
        res_instlist.append( [op, ili, oli] )

    return res_instlist, imm_instlist    
    
def fma_solver(instlist):
    regtype = lambda reg: reg.split("#")[0]
    
    #auxiliary register check
    aux_reg_list = []
    r_cnt, w_cnt = {},{}
    for line in instlist:
        op, reads, writes = line
        for reg in reads:
            r_cnt[reg] = r_cnt.setdefault(reg,0)+1
        for reg in writes:
            w_cnt[reg] = w_cnt.setdefault(reg,0)+1
    for reg, cnt in w_cnt.items():
        if regtype(reg) != "L":
            continue
        if cnt == 1 and r_cnt.setdefault(reg,0)<=1: # temporal register found
            aux_reg_list.append(reg)
    
    remain_lines = [ it for it in instlist ][::-1]
    ret_instlist = []
    fma_cnt = 0
    for i in range(10000):
        if len(remain_lines)==0:
            break
        line = remain_lines.pop()
        op, reads, writes = line
        if len(op.split(':')) < 2:
            ret_instlist.append(line)
            continue
        op_code, tp = op.split(':')
        if len(remain_lines)>0:
            opN, readsN, writesN = remain_lines[-1]
            if len(opN.split(':')) < 2:
                ret_instlist.append(line)
                continue
            opN_code, tpN = opN.split(':')
            if op_code == "mul" and opN_code in ["add","sub"] and tp == tpN:
                if writes[0] in aux_reg_list:
                    if writes[0] in readsN:
                        DC = [reg for reg in readsN if reg != writes[0] ][0]
                        DA, DB = reads
                        remain_lines.pop()
                        if opN_code == "add": # a*b + c or c + a*b
                            fma_op_name = "fma"
                        if opN_code == "sub" and writes[0] == readsN[0]: # a*b - c
                            fma_op_name = "fms"
                        if opN_code == "sub" and writes[0] == readsN[1]: # c - a*b
                            fma_op_name = "fmi"
                        ret_instlist.append(( fma_op_name+':'+tp ,[DA,DB,DC],writesN))
                        continue
        ret_instlist.append( line )
    return ret_instlist

def unvectorizeInst(instlist):
    rinst = []
    for op_with_type, ili, oli in instlist:
        if ":" in op_with_type:
            op, opt = op_with_type.split(":")
            opt = opt.split("vec")[0]
        else:
            op = op_with_type
            opt = ""
        #print (op, ili, "->", oli)
        if op == "mul":
            i1, i2 = ili
            o1, = oli
            i1n, i1t = i1.split(":")
            i2n, i2t = i2.split(":")
            o1n, o1t = o1.split(":")
            if ("vec" in i1t) and ("vec" in i2t): # dot product
                ut1 = i1t.split("vec")[0]
                ut2 = i2t.split("vec")[0]
                assert ut1 == ut2
                ut = ut1
                uto = o1t.split("vec")[0]
                rinst.append ([ f"mul:{opt}", [ i1n+"@x:"+ut, i2n+"@x:"+ut], [ o1n+"@a:"+ut ] ])
                rinst.append ([ f"fma:{opt}", [ i1n+"@y:"+ut, i2n+"@y:"+ut, o1n+"@a:"+ut], [ o1n+"@b:"+ut ] ])
                rinst.append ([ f"fma:{opt}", [ i1n+"@z:"+ut, i2n+"@z:"+ut, o1n+"@b:"+ut], [ o1n+":"+uto ] ])
                continue
            if ("vec" not in i1t) and ("vec" in i2t): # left broadcast elementwise multiple
                ut1 = i1t
                ut2 = i2t.split("vec")[0]
                ut = ut1
                uto = o1t.split("vec")[0]
                rinst.append ([ f"mul:{opt}", [ i1n+":"+ut1, i2n+"@x:"+ut2], [ o1n+"@x:"+uto ] ])  
                rinst.append ([ f"mul:{opt}", [ i1n+":"+ut1, i2n+"@y:"+ut2], [ o1n+"@y:"+uto ] ])  
                rinst.append ([ f"mul:{opt}", [ i1n+":"+ut1, i2n+"@z:"+ut2], [ o1n+"@z:"+uto ] ])
                continue
            if ("vec" in i1t) and ("vec" not in i2t): # right broadcast elementwise multiple
                ut1 = i1t.split("vec")[0]
                ut2 = i2t
                ut = ut2
                uto = o1t.split("vec")[0]
                rinst.append ([ f"mul:{opt}", [ i1n+"@x:"+ut1, i2n+":"+ut2], [ o1n+"@x:"+uto ] ])  
                rinst.append ([ f"mul:{opt}", [ i1n+"@y:"+ut1, i2n+":"+ut2], [ o1n+"@y:"+uto ] ])  
                rinst.append ([ f"mul:{opt}", [ i1n+"@z:"+ut1, i2n+":"+ut2], [ o1n+"@z:"+uto ] ]) 
                continue
        if op in ["fma","fms","fmi"]:
            i1, i2, i3 = ili # o1 = i1*i2 (+/-) i3
            o1, = oli
            i1n, i1t = i1.split(":")
            i2n, i2t = i2.split(":")
            i3n, i3t = i3.split(":")
            o1n, o1t = o1.split(":")
            if ("vec" in i1t) and ("vec" in i2t) and ("vec" in i3t) : # (vecA dot vecB)[scalar] + vecC <-- 特殊すぎる．．．
                assert "vec" in o1t # vector output required
                ut1 = i1t.split("vec")[0]
                ut2 = i2t.split("vec")[0]
                ut = ut1
                ut3 = i3t.split("vec")[0]
                uto = o1t.split("vec")[0]
                rinst.append ([ f"mul:{opt}", [ i1n+"@x:"+ut1, i2n+"@x:"+ut2], [ o1n+"@a:"+ut ] ])
                rinst.append ([ f"fma:{opt}", [ i1n+"@y:"+ut1, i2n+"@y:"+ut2, o1n+"@a:"+ut], [ o1n+"@b:"+ut ] ])
                rinst.append ([ f"fma:{opt}", [ i1n+"@z:"+ut1, i2n+"@z:"+ut2, o1n+"@b:"+ut], [ o1n+"@c:"+ut ] ])
                opX = {"fma":"add","fms":"sub"}[op]
                rinst.append ([ f"{opX}:{opt}", [ o1n+"@c:"+ut, i3n+"@x:"+ut3], [ o1n+"@x:"+uto ] ])
                rinst.append ([ f"{opX}:{opt}", [ o1n+"@c:"+ut, i3n+"@y:"+ut3], [ o1n+"@y:"+uto ] ])
                rinst.append ([ f"{opX}:{opt}", [ o1n+"@c:"+ut, i3n+"@z:"+ut3], [ o1n+"@z:"+uto ] ])
                continue
            if ("vec" in i1t) and ("vec" in i2t) and ("vec" not in i3t):
                assert "vec" not in o1t # scalar output required
                ut1 = i1t.split("vec")[0]
                ut2 = i2t.split("vec")[0]
                ut = ut1
                ut3 = i3t
                uto = o1t.split("vec")[0]
                ut = i3t       
                rinst.append ([ f"{op}:{opt}", [ i1n+"@x:"+ut1, i2n+"@x:"+ut2, i3n+":"+ut3], [ o1n+"@a:"+ut ] ])
                rinst.append ([ f"fma:{opt}",  [ i1n+"@y:"+ut1, i2n+"@y:"+ut2, o1n+"@a:"+ut], [ o1n+"@b:"+ut ] ])
                rinst.append ([ f"fma:{opt}",  [ i1n+"@z:"+ut1, i2n+"@z:"+ut2, o1n+"@b:"+ut], [ o1n+":"+uto ] ])
                continue
            if ("vec" not in i1t) and ("vec" in i2t) and ("vec" in i3t):
                assert "vec" in o1t
                ut1 = i1t
                ut2 = i2t.split("vec")[0]
                ut3 = i3t.split("vec")[0]
                uto = o1t.split("vec")[0]
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+"@x:"+ut2, i3n+"@x:"+ut3], [ o1n+"@x:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+"@y:"+ut2, i3n+"@y:"+ut3], [ o1n+"@y:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+"@z:"+ut2, i3n+"@z:"+ut3], [ o1n+"@z:"+uto ] ])
                continue
            if ("vec" in i1t) and ("vec" not in i2t) and ("vec" in i3t):
                assert "vec" in o1t
                ut1 = i1t.split("vec")[0]
                ut2 = i2t
                ut3 = i3t.split("vec")[0]
                uto = o1t.split("vec")[0]
                rinst.append ([ f"{op}:{opt}", [ i1n+"@x:"+ut1, i2n+":"+ut2, i3n+"@x:"+ut3], [ o1n+"@x:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+"@y:"+ut1, i2n+":"+ut2, i3n+"@y:"+ut3], [ o1n+"@y:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+"@z:"+ut1, i2n+":"+ut2, i3n+"@z:"+ut3], [ o1n+"@z:"+uto ] ])
                continue
            if ("vec" not in i1t) and ("vec" in i2t) and ("vec" not in i3t):
                assert "vec" in o1t
                ut1 = i1t
                ut2 = i2t.split("vec")[0]
                ut3 = i3t
                uto = o1t.split("vec")[0]
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+"@x:"+ut2, i3n+":"+ut3], [ o1n+"@x:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+"@y:"+ut2, i3n+":"+ut3], [ o1n+"@y:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+"@z:"+ut2, i3n+":"+ut3], [ o1n+"@z:"+uto ] ])
                continue
            if ("vec" in i1t) and ("vec" not in i2t) and ("vec" not in i3t):
                assert "vec" in o1t
                ut1 = i1t.split("vec")[0]
                ut2 = i2t
                ut3 = i3t
                uto = o1t
                rinst.append ([ f"{op}:{opt}", [ i1n+"@x:"+ut1, i2n+":"+ut2, i3n+":"+ut3], [ o1n+"@x:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+"@y:"+ut1, i2n+":"+ut2, i3n+":"+ut3], [ o1n+"@y:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+"@z:"+ut1, i2n+":"+ut2, i3n+":"+ut3], [ o1n+"@z:"+uto ] ])
                continue
            if ("vec" not in i1t) and ("vec" not in i2t) and ("vec" in i3t):
                assert "vec" in o1t
                ut1 = i1t
                ut2 = i2t
                ut3 = i3t.split("vec")[0]
                uto = o1t.split("vec")[0]
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+":"+ut2, i3n+"@x:"+ut3], [ o1n+"@x:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+":"+ut2, i3n+"@y:"+ut3], [ o1n+"@y:"+uto ] ])
                rinst.append ([ f"{op}:{opt}", [ i1n+":"+ut1, i2n+":"+ut2, i3n+"@z:"+ut3], [ o1n+"@z:"+uto ] ])
                continue
    
        if  len(oli) == 1: 
            o1n, o1t = oli[0].split(":")
            if ("vec" in o1t): #elementwise operator. (default vec-op unvectorizer)
                assert all(["vec" in il for il in ili])
                ut = o1t.split("vec")[0]
                rinst.append ([ f"{op}:{opt}", [ n+"@x:"+t.split("vec")[0] for n,t in [ it.split(":") for it in ili]  ] ,[ o1n+"@x:"+ut ]  ]) 
                rinst.append ([ f"{op}:{opt}", [ n+"@y:"+t.split("vec")[0] for n,t in [ it.split(":") for it in ili]  ] ,[ o1n+"@y:"+ut ]  ]) 
                rinst.append ([ f"{op}:{opt}", [ n+"@z:"+t.split("vec")[0] for n,t in [ it.split(":") for it in ili]  ] ,[ o1n+"@z:"+ut ]  ]) 
                continue
        if True: # others
            rinst.append ([ f"{op}:{opt}", ili, oli])
    return rinst

def generate_instlist_with_undetermined_raw_registers(core_instlist,decllist, imm_instlist):
    varlist_force= [ it[1][1] for decl_type, it in decllist if decl_type in ["decl_force", "decl_local_force"] ]
    varlist_epi= [ it[1][1] for decl_type, it in decllist if decl_type in ["decl_epi","decl_local_epi"] ]
    varlist_epj= [ it[1][1] for decl_type, it in decllist if decl_type in ["decl_epj","decl_local_epj"] ]
    
    #print ( "!f", varlist_force )
    #print ( "!i", varlist_epi )
    #print ( "!j", varlist_epj )
    
    
    rlist = []
    for op, ili, oli in core_instlist:
        rlist += ili
        rlist += oli
    rlist = set(rlist)
    
    rdic = {}
    roff = 0
    for i,it in enumerate(rlist):
        vartype = it.split(":")[1]
        if it.split("#")[0] == "M" :
            varname = it.split(":")[0].split("@")[0].split("#")[-1]
            name = it.split(":")[0].split("#")[-1]
            if varname in varlist_force:
                rdic[it] = f"F.{name}:{vartype}"
                continue
            if varname in varlist_epi:
                rdic[it] = f"I.{name}:{vartype}"
                continue
            if varname in varlist_epj:
                rdic[it] = f"J.{name}:{vartype}"
                continue
        if it.split("#")[0] in ["VF","VI",'C']:
            rdic[it] = it
            continue
        if it.split("#")[0] in ["O"]:
            rdic[it] = "O." + it.split('#')[1].split(':')[0] + f":{vartype}"
            continue
        if True: #Default
            rdic[it] = f"X.{roff}:{vartype}"
            roff +=1
            continue
        print ("?unsolved")
    #print (rdic)
   
    rinst = []
    for op_with_type, ili, oli in core_instlist:
        li_type = [ it.split(":")[-1] for it in ili ]
        op, tp = op_with_type.split(':')
        # F64 mul/fma is divided into two
        if op in ["mul","fma","fms","fmi"] and tp == "F64":
            rinst.append([ op+"u:"+tp, [ rdic[it] for it in ili ], [ rdic[it] for it in oli ] ])
            if op == "mul":
                rinst.append([ "fmad:"+tp, [ rdic[it] for it in ili ] + [ rdic[it] for it in oli ], [ rdic[it] for it in oli ] ])
            else:
                ab = [rdic[it] for it in ili][0:2]
                c = [ rdic[it] for it in oli ]
                next_op = "fmi" if op == "fmi" else "fma"
                rinst.append([ next_op+"d:"+tp, ab + c, [ rdic[it] for it in oli ] ])
        else:
            rinst.append([ op+":"+tp, [ rdic[it] for it in ili ], [ rdic[it] for it in oli ] ])

    #imm transfer:
    imm_instlist_gen = []
    for op, [val], [reg] in imm_instlist:
        #print ("imm", val, reg)
        imm_instlist_gen.append( [rdic[reg], val ] ) 
        
        
    return rinst, roff, len(rinst), imm_instlist_gen

def split_f64_for_f32(instlist):
    ret = []
    for inst in instlist:
        op,reads,writes = inst
        if op in ["mask","mend"]:
            ret.append(inst)
            continue
        opt = op.split(":")[1]
        wts = [ get_reg_type(reg) for reg in writes ]
        rts = [ get_reg_type(reg) for reg in reads ]
        assert all([w == wts[0] or w == "BOOL" for w in wts])
        isF64op = any([ r == "F64" for r in rts ])
        if isF64op and (wts[0] == "F32" or any([ r == "F32" for r in rts])):
            ret.append([op.replace("F32","F64"),[ r + ":0" for r in reads ], [ w + ":0" for w in writes]])
            ret.append([op.replace("F32","F64"),[ r + ":1" for r in reads ], [ w + ":1" for w in writes]])
        else:
            ret.append(inst)
    # check if mask with single/double mixed precision
    maskdepth = 0
    mask_tp = "F32"
    for inst in ret:
        if inst[0] == "mask:":
            maskdepth += 1
            continue
        if inst[0] == "mend:":
            maskdepth -= 1
            continue
        op, tp = inst[0].split(':')
        if op in ["cmp","cmpeq"]:
            mask_tp = tp
            continue
        if maskdepth > 0:
            if tp != mask_tp:
                print("error: multiprecision in if statement is not allowed")
                raise AssertionError()
    return ret

def get_reg_name(reg):
    return reg.split('.')[0]
def get_reg_id(reg):
    return reg.split(':')[0].split('.')[1]
def get_reg_type(reg):
    return reg.split(':')[1]

def add_rotate_for_simd(instlist,decllist,immlist,Rnum,Ilen):
    require_rotate = False
    simd_vars=[]
    regid=0
    for op, reads, writes in instlist:
        fp=op.split(':')[1]
        if fp == 'F32':
            for reg in reads:
                if 'J.' in reg:
                    if reg.split(':')[1] == "F32":
                        simd_vars.append((reg,fp))
                        require_rotate = True
        for reg in (reads+writes):
            if 'X.' in reg:
                reg = int(get_reg_id(reg))
                regid = max(reg,regid)
    regid = regid + 1
    rotate_insts=[]
    if require_rotate:
        Rnum = Rnum + 1
        for var, fp in simd_vars:
            if fp == 'F32':
                rotate_insts.append(['rotate:F32', [f'{var}', f'X.{regid}:F32'], [f'{var}']]) 
                Ilen = Ilen + 1
            if fp == 'F32vec':
                for dim in ['x','y','z']:
                    rotate_insts.append(['rotate:F32', [f'{var}',f'X.{regid}:F32'], [f'{var}']])
                Ilen = Ilen + 3
                requre_rotate = True
        immlist.append([f'X.{regid}:{fp}','VI#32:U64'])
    return rotate_insts+instlist, Rnum, Ilen, immlist
