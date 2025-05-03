def remove_redundant_set(ilist, decllist):
    order_graph().run(ilist, decllist)

####
class order_graph:
    def __init__(me, verbose=False):
        me.mask_stack = ['M']
        me.values = []
        me.graph = [] # node: (val_id:true, val_id:false, val_id:cond)
        me.value = 0
        me.verbose = verbose
    def get_value(me):
        me.value += 1
        return me.value
    def run(me, ilist, dlist):
        varlist_global = [ "M#"+it[1][1] for decl_type, it in dlist]
        varlist_global += [ ili[0] for op, ili, oli in ilist if op=="imm" ]
                 
        #Set pair Prepared
        pair_set = []
        pair_set_flag = []
        for pos, (op, ili, oli) in enumerate(ilist):
            if op == "set":
                pair_set.append((ili[0], oli[0]))
                pair_set_flag.append(True)
        #print (pair_set)
        
        #Set pair checking
        
        mask_stack = ["M"]
        RegVals  = {it:me.get_value() for it in varlist_global}
        RegVals_Start = RegVals.copy()
        MaskVals = {it:mask_stack.copy() for it in varlist_global}
        RegVals_hist = [] #[ RegVals.copy() ]
        for pos, (op, ili, oli) in enumerate(ilist):
            ri = ili[0] if len(ili)>0 else None
            ro = oli[0] if len(oli)>0 else None
            if op == "mask":
                mask_stack.append( ri )
            elif op == "mend":
                mask_stack.pop()
            elif op == "set":
                if RegVals.get(ro,None)==None or MaskVals[ri][-1] == MaskVals[ro][-1] or MaskVals[ri]==mask_stack:
                    #print ("&",op,ili,oli, MaskVals[ri][-1],mask_stack[-1])
                    RegVals[ro] = RegVals[ri]
                else:
                    RegVals[ro] = me.get_value()
                MaskVals[ro] = mask_stack.copy()
            else:
                RegVals[ro] = me.get_value()
                MaskVals[ro] = mask_stack.copy()
                
            RegVals_hist.append( RegVals.copy() )
        
        uniset = []
        for pos, (op, ili, oli) in enumerate(ilist):
            if op == "set":
                ri, ro = ili[0], oli[0]
                if RegVals_hist[pos][ri] == RegVals_hist[pos][ro]:
                    findset = [ k for k in uniset if (ri in k) or (ro in k) ]
                    if len(findset) == 0:
                        uniset.append( {ri,ro} )
                    else:
                        findset[0].add(ri)
                        findset[0].add(ro)
        remp = {}
        for fset in uniset:
            fli = sorted(list(fset))[::-1]
            #print (fli)
            remp.update({ it:fli[0] for it in fli[1:]})       
        print (remp) if me.verbose else None

        vlist = []
        for pos, (op, ili_, oli_) in enumerate(ilist):
            ili = [ remp.get(reg,reg) for reg in ili_] 
            oli = [ remp.get(reg,reg) for reg in oli_] 
            #print (op, ili, oli)
            if op=="set" and ili[0] == oli[0]:
                continue
            vlist.append([op, ili, oli])
        
        return vlist
