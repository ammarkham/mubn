
#include "io430.h"
#include "m_defs.h"
#include "m_arith.h"

void test_sum() {
    int errors = 0;
    uint16_t a_0[12] = {0x7DA0, 0x6C9D, 0x574F, 0xD615, 0x85C2, 0xCB5F, 0xCA49, 0x4285, 0x1F5C, 0x7799, 0xEA25, 0xD94C};
    uint16_t b_0[12] = {0xF3EA, 0x85D3, 0x0C32, 0xA075, 0xC431, 0x4AA2, 0xBC82, 0x96C5, 0xF7A7, 0x2B89, 0x0215, 0xD845};
    uint16_t c_0[12] = {0x718A, 0xF271, 0x6381, 0x768A, 0x49F4, 0x1602, 0x86CC, 0xD94B, 0x1703, 0xA323, 0xEC3A, 0xB191};
    uint16_t d_0[12];
    add_mp_elements(d_0, a_0, b_0, 12);
    if (0 == are_mp_equal(c_0, d_0, 12)) {
        errors++;
    }
}

int main( void )
{
  // Stop watchdog timer to prevent time out reset
  WDTCTL = WDTPW + WDTHOLD;
  test_sum();
  return 0;
}
