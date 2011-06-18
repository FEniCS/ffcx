class ArgumentNameFormatter(object):

    def name(self,argname):
        return argname

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
        'name',
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
        name = name or count or str(argument)

        # cache names of derived C++ types and objects
        self.names = dict((n,getattr(name_formatter,n)(name)) for n in self._names)

        # also save names as members for easy accessibility
        for k, v in self.names.iteritems():
            setattr(self,k,v)

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

    def evaluate_basis(self, pos):     
        code = '''
std::vector<{range_type}> {name}_basis({size_name});
{fe_switch_type}::basis({lfs_name}.finiteElement()).evaluateFunction({pos},{name}_basis);
'''
        names = self.names.copy()
        names['pos'] = pos
        return code.format(**names)
                    

    def evaluate_value(self, pos, iter_var='i'):
        from uflacs.codeutils.format_code import Block
        names = self.names.copy()
        names['iter_var'] = iter_var
        
        code = [
            '{range_type} {name}(0);'.format(**names),
            'for (std::size_t {iter_var} = 0; {iter_var} < {size_name}; ++{iter_var})'.format(**names),
            Block(
                '{name}.axpy(x({lfs_name},{iter_var}),{name}_basis[{iter_var}]);'.format(**names)
                )
            ]
        return code
