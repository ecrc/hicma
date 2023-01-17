/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file global.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global coreblas variables and functions
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @date 2010-11-15
 *
 */
static int coreblas_gemm3m_enabled = 0;

void
HICMA_set_coreblas_gemm3m_enabled( int v ) {
    coreblas_gemm3m_enabled = v;
}

int
HICMA_get_coreblas_gemm3m_enabled(void) {
    return coreblas_gemm3m_enabled;
}

/**
 *  LAPACK Constants
 */
char *hicma_lapack_constants[] =
        {
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "",                     // 100
                "Row",                  // 101: HicmaRowMajor
                "Column",               // 102: HicmaColMajor
                "", "", "", "", "", "", "", "",
                "No transpose",         // 111: HicmaNoTrans
                "Transpose",            // 112: HicmaTrans
                "Conjugate transpose",  // 113: HicmaConjTrans
                "", "", "", "", "", "", "",
                "Upper",                // 121: HicmaUpper
                "Lower",                // 122: HicmaLower
                "All",                  // 123: HicmaUpperLower
                "", "", "", "", "", "", "",
                "Non-unit",             // 131: HicmaNonUnit
                "Unit",                 // 132: HicmaUnit
                "", "", "", "", "", "", "", "",
                "Left",                 // 141: HicmaLeft
                "Right",                // 142: HicmaRight
                "", "", "", "", "", "", "", "",
                "",                     // 151:
                "",                     // 152:
                "",                     // 153:
                "",                     // 154:
                "",                     // 155:
                "",                     // 156:
                "Epsilon",              // 157: HicmaEps
                "",                     // 158:
                "",                     // 159:
                "",                     // 160:
                "", "", "", "", "", "", "", "", "", "",
                "One norm",             // 171: HicmaOneNorm
                "",                     // 172: HicmaRealOneNorm
                "",                     // 173: HicmaTwoNorm
                "Frobenius norm",       // 174: HicmaFrobeniusNorm
                "Infinity norm",        // 175: HicmaInfNorm
                "",                     // 176: HicmaRealInfNorm
                "Maximum norm",         // 177: HicmaMaxNorm
                "",                     // 178: HicmaRealMaxNorm
                "",                     // 179
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "",                     // 200
                "Uniform",              // 201: HicmaDistUniform
                "Symmetric",            // 202: HicmaDistSymmetric
                "Normal",               // 203: HicmaDistNormal
                "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "",                     // 240
                "Hermitian",            // 241 HicmaHermGeev
                "Positive ev Hermitian",// 242 HicmaHermPoev
                "NonSymmetric pos sv",  // 243 HicmaNonsymPosv
                "Symmetric pos sv",     // 244 HicmaSymPosv
                "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "",                     // 290
                "No Packing",           // 291 HicmaNoPacking
                "U zero out subdiag",   // 292 HicmaPackSubdiag
                "L zero out superdiag", // 293 HicmaPackSupdiag
                "C",                    // 294 HicmaPackColumn
                "R",                    // 295 HicmaPackRow
                "B",                    // 296 HicmaPackLowerBand
                "Q",                    // 297 HicmaPackUpeprBand
                "Z",                    // 298 HicmaPackAll
                "",                     // 299

                "",                     // 300
                "No vectors",           // 301 HicmaNoVec
                "Vectors needed",       // 302 HicmaVec
                "I",                    // 303 HicmaIvec
                "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "", "", "", "", "", "", "", "", "", "",
                "",                     // 390
                "Forward",              // 391
                "Backward",             // 392
                "", "", "", "", "", "", "", "",
                "Columnwise",           // 401
                "Rowwise",              // 402
                "", "", "", "", "", "", "", ""  // Remember to add a coma!
        };
