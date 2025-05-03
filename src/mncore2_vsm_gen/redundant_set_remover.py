def redundant_set_remover(instlist, decllist):
    #instlist, decllist = tree_trace_result.instlist ,tree_trace_result.decllist
    instlist = reduceEquivalentReg(instlist, decllist)
    #print ("Round 2-----------------")
    instlist = reduceEquivalentReg(instlist, decllist)
    return instlist

    
def reduceEquivalentReg(ilist,decllist):
    varlist_global = [ "M#"+it[1][1] for decl_type, it in decllist]
    #print( "###GLOBAL VAR###:", varlist_global )
    
    def subst_simulate( nlist, regA, regB ):
        #print ("****", regA, regB)
        
        RegVal = { it: it for v,it in enumerate(varlist_global) }
        mask_stack = []
        
        for pos, (op, ili, oli) in enumerate(nlist):
            #print (op,ili,"->",oli)

            for ri in ili:
                ri = regB if ri == regA else ri
                if ri not in RegVal:
                    RegVal[ri] = ri
                    

            def masked( mask_stack, valnew, valold ):
                maskstr = ",".join( mask_stack )
                if valold == "_" or valold == valnew :
                    return valnew
                return f'[{maskstr}?{valnew}:{valold}]'
                
            if op == "mask":
                ri = ili[0]
                ri = regB if ri == regA else ri
                mask_stack.append( RegVal[ri] )
            elif op == "mend":
                mask_stack.pop()

            elif op == "set":
                ri, ro = ili[0], oli[0]
                ri = regB if ri == regA else ri
                ro = regB if ro == regA else ro
                RegVal[ro] = masked( mask_stack, RegVal[ri] ,RegVal.get(ro,"_") )
                
            elif op == "imm":
                ri, ro = ili[0], oli[0]
                ri = regB if ri == regA else ri
                ro = regB if ro == regA else ro
                if ro[:2] == "L#":
                    RegVal[ro] = masked( mask_stack, RegVal[ri] ,RegVal.get(ro,"_") )
                else:
                    RegVal[ro] = "ERR!"
                    #print ("ERR!")

            elif len(oli)>0:
                ro = oli[0]
                ro = regB if ro == regA else ro
                val = op +"(" + ",".join( [RegVal[regB if ri == regA else ri] for ri in ili] ) +")"
                RegVal[ro] = masked( mask_stack, val ,RegVal.get(ro,"_") )
                #print ("\t\t", ro, RegVal[ro] ) if regA == None else None
            #print ( "\t", op,ili,"->",oli, RegVal ) if regA == None else None
                        
        if regA != regB and regA is not None and regB is not None:
            if regA in RegVal and regB in RegVal:
                RegVal[regA] = RegVal[regB]
            
        return { key:RegVal[key] for key in varlist_global }
                        

    RegVal_Main = subst_simulate( ilist,None,None )
    #print ("MAIN")
    #print (RegVal_Main)
        
    rlist = [ (op, ili, oli) for op, ili, oli in ilist ]
    for pos, (op, ili, oli) in enumerate(ilist):
        if op == "set":
            ri, ro = sorted([ili[0], oli[0]],key=lambda reg: (3 if reg in varlist_global else 0) + (2 if reg[:2] == "M#" else 0) + (1 if reg[0] == "F" else 0)  )
            RegVal_sub = subst_simulate( rlist, ri, ro )
            #print ("?")
            #print (RegVal_sub)
            if RegVal_sub == RegVal_Main:
                #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@set removed", ri, ro)
                rlist = [ (op, [ ro if reg==ri else reg for reg in ili] , [ ro if reg==ri else reg for reg in oli]) for op, ili, oli in rlist ]
                #print (RegVal_sub)
                #print_instlist(rlist)
                continue
                
            ri, ro = ro, ri
            RegVal_sub = subst_simulate( rlist, ri, ro )
            #print ("?")
            #print (RegVal_sub)
            if RegVal_sub == RegVal_Main:
                #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@set removed", ri, ro)
                rlist = [ (op, [ ro if reg==ri else reg for reg in ili] , [ ro if reg==ri else reg for reg in oli]) for op, ili, oli in rlist ]
                #print (RegVal_sub)
                #print_instlist(rlist)
                continue
            else:
                #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@not not not", ri, ro)
                #print (RegVal_sub)
                #print_instlist(rlist)
                continue
    
    #Remove loop
    tlist = []
    for pos, (op, ili, oli) in enumerate(rlist):
        if op == "set" and ili[0]==oli[0]:
            continue
        tlist.append( [op, ili, oli] )
    
    return tlist