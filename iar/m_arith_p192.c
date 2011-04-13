
/*
 * Copyright 2010 UMass Amherst. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice, this list of
 *       conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list
 *       of conditions and the following disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY UMass Amherst ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
//#include <__cross_studio_io.h>
//#include "msp430setup.h"
#include <signal.h>
#include <stdio.h>
#include <string.h>

#include "m_defs.h"
#include "m_arith.h"
#include "m_arith_p192.h"

#ifdef DEBUG_PRINTF
#include "debug_utils.h"
#endif
/**
 * Fast reduction modulo p = p_192 = 2^{192} - 2^{64} -1 c.f Alg. 2.27
 * Input: c, s.t 0<= c <= p^2
 * Output: c mod p
 */
void reduce_mod_p(uint16_t * c, uint16_t * p, uint16_t * c_out) {
    uint16_t s1[12];
    uint16_t s2[12];
    uint16_t s3[12];
    uint16_t s4[12];

    s1[0] = c[0 + 0];
    s1[1] = c[1 + 0];
    s1[2] = c[2 + 0];
    s1[3] = c[3 + 0];
    s1[4] = c[0 + 4];
    s1[5] = c[1 + 4];
    s1[6] = c[2 + 4];
    s1[7] = c[3 + 4];
    s1[8] = c[0 + 8];
    s1[9] = c[1 + 8];
    s1[10] = c[2 + 8];
    s1[11] = c[3 + 8];

    s2[0] = c[0 + 12];
    s2[1] = c[1 + 12];
    s2[2] = c[2 + 12];
    s2[3] = c[3 + 12];
    s2[4] = c[0 + 12];
    s2[5] = c[1 + 12];
    s2[6] = c[2 + 12];
    s2[7] = c[3 + 12];
    s2[8] = 0;
    s2[9] = 0;
    s2[10] = 0;
    s2[11] = 0;

    s3[0] = 0;
    s3[1] = 0;
    s3[2] = 0;
    s3[3] = 0;
    s3[4] = c[0 + 16];
    s3[5] = c[1 + 16];
    s3[6] = c[2 + 16];
    s3[7] = c[3 + 16];
    s3[8] = c[0 + 16];
    s3[9] = c[1 + 16];
    s3[10] = c[2 + 16];
    s3[11] = c[3 + 16];

    s4[0] = c[0 + 20];
    s4[1] = c[1 + 20];
    s4[2] = c[2 + 20];
    s4[3] = c[3 + 20];
    s4[4] = c[0 + 20];
    s4[5] = c[1 + 20];
    s4[6] = c[2 + 20];
    s4[7] = c[3 + 20];
    s4[8] = c[0 + 20];
    s4[9] = c[1 + 20];
    s4[10] = c[2 + 20];
    s4[11] = c[3 + 20];

    //return s1 + s2 + s3 +s4

    add_mod_p(c_out, s3, s4, p, 12);
    add_mod_p(c_out, c_out, s2, p, 12);
    add_mod_p(c_out, c_out, s1, p, 12);
}

void multiply_mod_p_192(uint16_t * c, uint16_t * a, uint16_t * b) {
    /*
     * Using NIST prime p_192 = 2^{192} - 2^{64} -1
     */
    uint16_t p[WORDLENGTH] = {0xffff, 0xffff, 0xffff, 0xffff, 0xfffe, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff, 0xffff};

    uint16_t out[2*WORDLENGTH];

    multiply_mp_elements(out, a, b, 12);
    //print_bn((uint8_t *) "c", out, 24);
    reduce_mod_p(out, p, c);
    //print_bn((uint8_t *) "c mod p", out, 12);
}

/**
 * Implementation of the left to right modular exponentiation algorithm
 * as described in the HAC book by Menezes et.al.
 *
 * @param A The result of raising g to the power of e
 * @param g an element of Z*_p
 * @param e a single precission exponent
 */
void mod_exp_p_192_lr(uint16_t * A, uint16_t * g, uint16_t e) {
    uint16_t temp[12];
    int i;
    int t = bit_length(e);

    //1.
#ifdef DEBUG_PRINTF
    printf("1.\n");
#endif
    set_to_zero(A, 12);
    A[0] = 1;
#ifdef DEBUG_PRINTF
    print_bn((uint8_t *) "A", A, 12);
    print_bn((uint8_t *) "g", g, 12);
#endif
    //2.

    for (i = t; i >= 0; i--) { // Note, first decrease, then work
#ifdef DEBUG_PRINTF
        printf("exp i: %d\n", i);
#endif
        //2.1 A = A*A mod p
        multiply_mod_p_192(temp, A, A);
#ifdef DEBUG_PRINTF
        print_bn((uint8_t *) "A*A mod p", temp, 12);
#endif
        copy_mp(A, temp, 12);
        //2.2 If e_i = 1 then A = Mont(A,x_hat)
        if (1 == ith_bit(e, i)) {
#ifdef DEBUG_PRINTF
            printf("%d-bit = 1\n", i);
#endif
            multiply_mod_p_192(temp, A, g);
#ifdef DEBUG_PRINTF
            print_bn((uint8_t *) "A*g", temp, 12);
#endif
            copy_mp(A, temp, 12);
#ifdef DEBUG_PRINTF
            printf("exting here?\n");
#endif
        }
#ifdef DEBUG_PRINTF
        printf("or exting here?\n");
#endif
    }
    /*
    for (i = t - 1; i >= 0; i = i-1){
        printf("exp i: %d\n", i);
        //2.1 A = A*A mod p
        multiply_mod_p_192(temp, A, A);
        print_bn((uint8_t *) "A*A mod p", temp, 12);
        copy_mp(A, temp, 12);
    }*/

    //3.
#ifdef DEBUG_PRINTF
    printf("result:\n");
#endif
    //copy_mp(out,A,12);

}

/**
 * Implementation of the left to right modular exponentiation algorithm
 * as described in the HAC book by Menezes et.al.
 *
 * @param A The result of raising g to the power of e
 * @param g an element of Z*_p
 * @param e a multi-precission exponent
 * @param e_legth the wordlength of the multi-precission exponent
 */

void mod_exp_p_192(uint16_t * A, uint16_t * g, uint16_t * e, uint16_t e_length) {
    uint16_t temp[12];
    int i;
    int t = mp_bit_length(e,e_length);

    //1.
#ifdef DEBUG_PRINTF
    printf("1.\n");
#endif
    set_to_zero(A, 12);
    A[0] = 1;
#ifdef DEBUG_PRINTF
    print_bn((uint8_t *) "A", A, 12);
    print_bn((uint8_t *) "g", g, 12);
#endif
    //2.

    for (i = t; i >= 0; i--) { // Note, first decrease, then work
#ifdef DEBUG_PRINTF
        printf("exp i: %d\n", i);
#endif
        //2.1 A = A*A mod p
        multiply_mod_p_192(temp, A, A);
#ifdef DEBUG_PRINTF
        print_bn((uint8_t *) "A*A mod p", temp, 12);
#endif
        copy_mp(A, temp, 12);
        //2.2 If e_i = 1 then A = Mont(A,x_hat)
        if (1 == mp_ith_bit(e, i)) {
#ifdef DEBUG_PRINTF
            printf("%d-bit = 1\n", i);
#endif
            multiply_mod_p_192(temp, A, g);
#ifdef DEBUG_PRINTF
            print_bn((uint8_t *) "A*g", temp, 12);
#endif
            copy_mp(A, temp, 12);
#ifdef DEBUG_PRINTF
            printf("exting here?\n");
#endif
        }
#ifdef DEBUG_PRINTF
        printf("or exting here?\n");
#endif
    }
    /*
    for (i = t - 1; i >= 0; i = i-1){
        printf("exp i: %d\n", i);
        //2.1 A = A*A mod p
        multiply_mod_p_192(temp, A, A);
        print_bn((uint8_t *) "A*A mod p", temp, 12);
        copy_mp(A, temp, 12);
    }*/

    //3.
#ifdef DEBUG_PRINTF
    printf("result:\n");
#endif
    //copy_mp(out,A,12);

}

