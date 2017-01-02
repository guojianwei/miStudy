# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 10:04:51 2016

@author: G2
"""
"""
INN-NB 模型和参数参考文献：
individual nearest-neighbor
Thermodynamic Parameters for an Expanded Nearest-Neighbor Model for Formation of RNA Duplexes with Watson−Crick Base Pairs
Figure 2. Calculation of thermodynamic properties for a non-self-complementary duplex,
:)
"""
class INN_HB:
    def __init__(self):
        self._delg_init = 4.09
        self._delg_tau = 	0.45
        self._delg_sym = 0
        bases1 = ['AA','AU','UA','CU','CA','GU','GA','CG','GG','GC']
        bases2 = ['UU','AU','UA','AG','UG','AC','UC','CG','CC','GC']
        dNN = [-0.93,-1.10,-1.33,-2.08,-2.11,-2.24,-2.35,-2.36,-3.26,-3.42]
        self._dict_NN = {}
        for i,v in enumerate(bases1):
            self._dict_NN[v] = dNN[i]
        for i,v in enumerate(bases2):
            self._dict_NN[v] = dNN[i]
            
    def base_pairing(self,bases53):
        '''
        < 从5'到3'的单链序列
        >inn-nb 的模型自由能
        '''
        bases53 = bases53.upper()
        sigma = 0
        for i in range(len(bases53)-1):
            mer2 = bases53[i:i+2]
            sigma += self._dict_NN[mer2]
            #print mer2,dict_NN[mer2]
        m_tau = 0
        if bases53[0] == 'A' or bases53[0] == 'U':
            m_tau += 1
        if bases53[-1] == 'A' or bases53[-1] == 'U':
            m_tau += 1
        delg_duplex = self._delg_init+ sigma + \
                        m_tau*self._delg_tau +self._delg_sym
        #print bases53,delg_duplex
        return delg_duplex

if __name__ == '__main__':
    bases53 = 'ACGCA'
    inn = INN_HB()
    print inn.base_pairing(bases53)