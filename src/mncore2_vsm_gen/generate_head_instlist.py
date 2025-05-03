from generate_core_instlist import get_reg_id, get_reg_name, get_reg_type

def split_reg(reg):
    items = reg.split(':')
    tp = items[1]
    simd = '' if len(items) < 3 else items[2]
    reg = items[0]
    splitted = reg.split(".")
    if len(splitted) == 2:
        regtype, regname = splitted
    elif len(splitted) == 3:
        regtype, iotype, regname = splitted
        regname = iotype + "." + regname
    return [regtype, regname, tp, simd]

def generate_head_instlist(core_instlist_al, state_al, imm_instlist_gen, decllist):
    regs_dic = gen_regs_dic(core_instlist_al)                                     
    core_instlist_head = gen_core_instlist_head(core_instlist_al, regs_dic)
    imm_instlist_head = gen_imm_instlist_head(imm_instlist_gen, regs_dic, state_al)
    
    save_obj = {
        "regs_dic": regs_dic,
        "core_instlist_head":core_instlist_head,
        "imm_instlist_head":imm_instlist_head,
        "decllist":decllist
    }
    return save_obj

##################
forwardings=["A","U","Z"]
isImm = lambda x: x.split("#")[0] in ["VF","VI","C"]
def gen_regs_dic(allocated_instruction):
    #solve registers
    regs_dic = { it:[] for it in "IJFRSMNO"}
    for insts in allocated_instruction:
        for op, reads, writes, md in insts:
            regs = reads + writes
            for reg in regs:
                if isImm(reg):
                    continue
                regtype,regname,tp,_ = split_reg(reg)
                if regtype in forwardings:
                    continue
                if f'{regname}:{tp}' not in regs_dic[regtype]:
                    regs_dic[regtype].append(f"{regname}:{tp}")
    #[ v.sort() for k,v in regs_dic.items() ]
    #[ print (k,v) for k,v in regs_dic.items() ]
    return regs_dic


def gen_core_instlist_head(allocated_instruction, regs_dic):
    instlist = []
    for insts in allocated_instruction:
        inst=[]
        for op, reads, writes, md in insts:
            reads_new = []
            for reg in reads:
                if isImm(reg):
                    reads_new.append ( reg )
                else:
                    regtype,regname,tp,simd = split_reg(reg)
                    if simd != "":
                        simd = ':' + simd 
                    if regtype in forwardings:
                        reads_new.append ( reg )
                    else:
                        reads_new.append ( f"{regtype}.{regs_dic[regtype].index(f'{regname}:{tp}')}:{tp}{simd}" )
            writes_new = []
            for reg in writes:
                regtype,regname, tp, simd = split_reg(reg)
                if simd != "":
                    simd = ':' + simd 
                writes_new.append ( f"{regtype}.{regs_dic[regtype].index(f'{regname}:{tp}')}:{tp}{simd}" if reg.split('.')[0] not in ['A','Z','U','O'] else reg )
            inst.append( [op, reads_new, writes_new, md] )
        instlist.append(inst)
    return instlist

def gen_imm_instlist_head(imm_instlist_gen, regs_dic, state_al):
    instlist = []
    rmn_li, t_li = state_al
    #print ( imm_instlist_gen, state_al )
    for reg, val in imm_instlist_gen:
        regname = split_reg(reg)[1]
        v = val.split("#")[1]
        loc = "rsmn"[rmn_li[int(regname)]]
        regtype = loc.upper()
        tp = get_reg_type(reg)
        instlist.append( [f"{regtype}.{regs_dic[regtype].index(f'{regname}:{tp}')}:{tp}", v ] )
    return instlist
