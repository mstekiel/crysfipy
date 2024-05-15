import crysfipy

if __name__=='__main__':
    cfp = crysfipy.CEFpars(dict(B21=1.203, B40=-0.001, B53=0.244), 'meV')
    print(cfp)
    cefion = crysfipy.CEFion(ion=crysfipy.Ion('Ce'),cfp=cfp, Hfield=[0,0,0])
    print(cefion)
