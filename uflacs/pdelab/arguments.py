class ArgumentNameFormatter(object):

    def lfs_type(self,argname):
        return 'LFS' + argname.upper()

    def lfs_name(self,argname):
        return 'lfs' + argname

    def fe_switch_type(self,argname):
        return 'FESwitch_' + self.lfs_type(argname)

    def basis_switch_type(self,argname):
        return 'BasisSwitch_' + self.lfs_type(argname)

    def size_name(self,argname):
        return self.lfs_name(argname) + '_size'

    def domain_field_type(self,argname):
        return 'DF_' + argname.upper()

    def range_field_type(self,argname):
        return 'RF_' + argname.upper()

    def range_type(self,argname):
        return 'R_' + argname.upper()


class ArgumentData(object):

    _names = [
        'lfs_type',
        'lfs_name',
        'fe_switch_type',
        'basis_switch_type',
        'size_name',
        'domain_field_type',
        'range_field_type',
        'range_type',
        ]

    def __init__(self,argument,name=None,count=None,name_formatter=ArgumentNameFormatter()):
        self.argument = argument
        self.name = name or count or str(argument)

        # cache names of derived C++ types and objects
        for n in self._names:
            setattr(self,n,getattr(name_formatter,n)(self.name))

        from uflacs.codeutils.format_code import TypeDef, Type

        # code for argument:
        code = []
        
        # typedef for FE interface switch
        code.append(
            TypeDef(
                Type('Dune::FiniteElementInterfaceSwitch',
                     ('typename %s::Traits::FiniteElementType' % self.lfs_type,),True
                     ),
                self.fe_switch_type)
            )
        
        # typedef for basis interface switch
        code.append(
            TypeDef(
                Type('Dune::BasisInterfaceSwitch',
                     ('typename %s::Basis' % self.fe_switch_type,),True
                     ),
                self.basis_switch_type)
             )   

        # typedef for domain field
        code.append(
            TypeDef(
                Type('typename %s::DomainField' % self.basis_switch_type),
                self.domain_field_type)
             )   
            
        # typedef for range field
        code.append(
            TypeDef(
                Type('typename %s::RangeField' % self.basis_switch_type),
                self.range_field_type)
             )   

        # typedef for range
        code.append(
            TypeDef(
                Type('typename %s::Range' % self.basis_switch_type),
                self.range_type)
             )

        # const for space size
        code.append(
            'const std::size_t %s = %s.size();' % (self.size_name, self.lfs_name)
            )

        self.preamble_code = code
