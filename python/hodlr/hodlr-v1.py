import sys
import numpy as np
from numpy import linalg as LA
from Tree import Tree
## WARNINGS: int() type casts must be removed
np.set_printoptions(suppress=False,precision=3, threshold=np.inf,  linewidth=10000)
__m__=16 #number of rows of dense _A_
__n__=__m__ #number of cols of dense _A_
__k__=1  #number of levels
__p__=2  #rank/number of columns of U

print("Script name: ", sys.argv[0])
print("Number of arguments: ", len(sys.argv))
print("Arguments: " , str(sys.argv))
if (len(sys.argv) == 4):
    __m__=int(sys.argv[1])
    __k__=int(sys.argv[2])
    __p__=int(sys.argv[3])
print("Usage: ")
print("\tm: number of rows of dense A")
#__n__=__m__ #number of cols of dense A
print("\tk: number of levels")
print("\tp: max rank/number of columns of U")


def debugprint(*args, print_option=True):
    if print_option is True:
        print(" ".join(map(str,args)))
def global_row_index(level, localindex, total_num_tiles, total_num_levels):
   gi=localindex*int(pow(2, total_num_levels-level))
   if gi >= total_num_tiles:
       return None
   return gi
def global_col_index(level, localindex, total_num_tiles, total_num_levels):
   gi=localindex*int(pow(2, total_num_levels-level))
   if gi >= total_num_tiles:
       return None
   return gi
def numtiles_level(level):
    return int(np.ceil(pow(2, level)))

print_init=False
print_init_diag=not True

_I_=np.eye(__m__)
_u_=np.random.random((__m__, __p__))
debugprint("u:\n", _u_, print_option=print_init or print_init_diag)

__t__=int(np.ceil(pow(2, __k__)))  #number of tiles
__D_m__=int(__m__/__t__)  #number of rows/cols of diagonal blocks
__U_m__=int(__m__/__t__)  #number of rows of U
__VT_n__=int(__n__/__t__)  #number of rows of VT

_D_=np.zeros(shape=(__m__ * __D_m__))
_U_=np.zeros(shape=(__m__ * __k__ *__p__))
_VT_=np.zeros(shape=(__m__ * __k__ *__p__))
ranks=np.zeros(shape=(__t__ * __t__), dtype=int) # OPTIMIZE: size of this array might be smaller
tree=Tree(__m__, __p__)
##################### INIT DIAGONALS
dststart=0
dstsize=__D_m__*__D_m__
drowstart=0
for t in range(0, __t__):
    debugprint("Tile:", t, print_option=print_init_diag)
    dstend  = dststart + dstsize
    drowend = drowstart + __D_m__
    _I2_ = _I_[drowstart:drowend, drowstart:drowend]
    debugprint("I:", _I2_, print_option=print_init_diag)
    _u2_ = _u_[drowstart:drowend]
    debugprint("u:", _u2_, print_option=print_init_diag)
    _u2_u2T_I_ = np.add(np.matmul(_u2_, _u2_.T), _I_[drowstart:drowend, drowstart:drowend])
    debugprint("uuT+I:\n", _u2_u2T_I_, print_option=print_init_diag)
    _D_[dststart:dstend] = _u2_u2T_I_.flatten('F')
    debugprint(_D_[dststart:dstend], print_option=print_init_diag)
    dststart  += dstsize
    drowstart += __D_m__
################### INIT OFF-DIAGONALS
for level in range(1, __k__+1):
    debugprint("Level:", level, print_option=print_init)
    tiles = Tree.tiles_at_depth(level)
    for tile in tiles:
        ranks[tile] = int(__p__)
        ststart  = tree.ststart(tile)
        stend    = tree.stend(tile)
        rowstart = tree.rowstart(tile)
        rowend   = tree.rowend(tile)
        colstart = tree.colstart(tile)
        colend   = tree.colend(tile)
        debugprint(tile, "  storage:", ststart, stend, "  row:", rowstart, rowend, "  col:", colstart, colend, print_option=print_init)
        _u2_ = _u_[rowstart:rowend]
        debugprint(_u2_, print_option=print_init)
        _u3_ = _u2_.flatten('F')
        _vt2_ = _u_[colstart:colend]
        debugprint(_vt2_, print_option=print_init)
        _vt3_ = _vt2_.flatten('F')
        _U_[ststart:stend] = _u3_
        _VT_[ststart:stend] =  _vt3_
        debugprint(_U_[ststart:stend], print_option=print_init)
        debugprint(_VT_[ststart:stend], print_option=print_init)
################## DECOMPRESS HODLR
print_decompress=True
print_decompress_diag=True
_dcAt_ = [[np.zeros(shape=(__D_m__, __D_m__)) for j in range(0, __t__)] for i in range(0, __t__)]
dststart=0
dstsize=__D_m__*__D_m__ #ASSUMPTION: square tiles
for t in range(0, __t__):
    debugprint("Tile:", t, print_option=print_decompress_diag)
    dstend  = dststart + dstsize
    debugprint(_D_[dststart:dstend], print_option=print_decompress_diag)
    _dcAt_[t][t] = np.reshape(_D_[dststart:dstend], (__D_m__, __D_m__))
    dststart  += dstsize
#TODO decompress dcAt into dcA
print("decompressed tile A:")
for i in range(0, __t__):
    for j in range(0, __t__):
         print(i,j)
         print(_dcAt_[i][j])
sys.exit()
for level in range(1, __k__+1):
    debugprint("Level:", level, print_option=print_decompress)
    tiles = Tree.tiles_at_depth(level)
    for tile in tiles:
        ststart  = tree.ststart(tile)
        stend    = tree.stend(tile)
        rowstart = tree.rowstart(tile)
        rowend   = tree.rowend(tile)
        colstart = tree.colstart(tile)
        colend   = tree.colend(tile)
        debugprint(tile, "  storage:", ststart, stend, "  row:", rowstart, rowend, "  col:", colstart, colend, print_option=print_decompress)
        # _u2_ = _u_[rowstart:rowend]
        #debugprint(_u2_, print_option=print_decompress)
        #_u3_ = _u2_.flatten('F')
        #_vt2_ = _u_[colstart:colend]
        #debugprint(_vt2_, print_option=print_decompress)
        #_vt3_ = _vt2_.flatten('F')
        _u2_ = np.reshape(_U_[ststart:stend], (rowend-rowstart, ranks[tile]))
        _vt2_ = np.reshape(_VT_[ststart:stend], (rowend-rowstart, ranks[tile]))
        print(_u2_)
        #_VT_[ststart:stend] =  _vt3_
        #debugprint(_U_[ststart:stend], print_option=print_decompress)
        #debugprint(_VT_[ststart:stend], print_option=print_decompress)
print("decompressed tile A:")
for i in range(0, __t__):
    for j in range(0, __t__):
         print(i,j)
         print(_dcAt_[i][j])
#sys.exit()
#norm_dcA2=LA.norm(_dcA2_, 2)
#print('norm(dcA)=', norm_dcA2)

print_dense_steps=not True
print_dense_A=True
_I_=np.eye(__m__)
debugprint("I:\n", _I_, print_option=print_dense_steps)
_u_uT_=np.matmul(_u_, _u_.T)
_A_ = np.add(_I_,_u_uT_)
debugprint("A:\n", _A_, print_option=print_dense_steps or print_dense_A)
norm_A=LA.norm(_A_, 2)
print('norm(A)=', norm_A)
sys.exit()
_L_=np.linalg.cholesky(_A_)
debugprint("L:\n", _L_, print_option=print_dense_steps)
_L_LT_=np.matmul(_L_, _L_.T)
debugprint("L*L^T:\n", _L_LT_, print_option=print_dense_steps)
_L_LT_minus_A_ = np.subtract(_L_LT_, _A_)
norm_diff=LA.norm(_L_LT_minus_A_, 2)
print('norm(A)=', norm_A, '  norm(LLT-A)=', norm_diff, '  norm(LLT-A)/norm(A)=', norm_diff/norm_A)

#https://docs.scipy.org/doc/scipy/reference/linalg.lapack.html
#https://docs.scipy.org/doc/scipy/reference/linalg.blas.html
sys.exit()
for level in range(1, __k__+1):
    numtiles=numtiles_level(level)
    numrows=int(__m__/numtiles)
    print("numtiles=", numtiles, "  having ", numrows, "rows    at level=", level)
    for li in range(0, numtiles):
        gi = global_row_index(level, li, __t__, __k__)
        usi = gi*__U_m__
        uei = gi*__U_m__+numrows
        gci = global_col_index(level, li, __t__, __k__)
        vtsi = 0
        vtei = 0
        print("level=", level, "  li=", li, "   gi=", gi, "  gci=", gci, "  u si-ei=", usi, uei, "  vt si-ei=", vtsi, vtei)
        _U_[level][gi] = _u_[usi:uei]
        _VT_[level][gi] = _u_[vtsi:vtei] ##ASSUMPTION __U_m__ = __VT_n__
        print(_U_[level][gi])
        print(_VT_[level][gi])
