import crysfipy

if __name__ == '__main__':
    # print(_implemented_StevensOps)
    Bdict = dict(B20=1, B52=0.1, B6m1=-0.3)
    cfp = crysfipy.CEFpars(Bdict, units='meV', pointGroup='2')

    print(cfp)

    print(crysfipy.CEFpars.allowed_Bpars(pointGroupName='4mm'))