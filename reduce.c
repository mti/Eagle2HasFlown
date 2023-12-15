/*
 * Implementors: EagleSign Team
 * This implementation is highly inspired from Dilithium and
 * Falcon Signatures' implementations
 */

#include <stdint.h>
#include "params.h"
#include "reduce.h"

/*
 * For finite field element y with -Q <= y < Q, compute r \equiv y(mod Q)
 * such that r is between -q/2 and +q/2.
 */
S_DOUBLE_Q_SIZE reduce(S_DOUBLE_Q_SIZE y)
{
  DOUBLE_Q_SIZE r;
  r = (DOUBLE_Q_SIZE)y;
  r += ((-Q) & ~-((r - ((Q >> 1) + 1)) >> (DOUBLE_Q_BIT_SIZE - 1))) | (Q & -((r + (Q >> 1)) >> (DOUBLE_Q_BIT_SIZE - 1)));
  return r;
}

/*
 * For small integer r between -q/2 and +q/2, compute the finite field
 * element y \equiv r(mod Q) with 0 <= y < Q
 */
Q_SIZE addq(S_Q_SIZE r)
{
  Q_SIZE y;
  y = (Q_SIZE)r;
  y += Q & -(y >> (Q_BIT_SIZE - 1));
  return y;
}

/*
 * Addition modulo q. Operands must be in the 0..q-1 range.
 */
S_DOUBLE_Q_SIZE add(S_DOUBLE_Q_SIZE x, S_DOUBLE_Q_SIZE y)
{
  DOUBLE_Q_SIZE d;

  d = x + y;
  d += ((-Q) & ~-((d - ((Q >> 1) + 1)) >> (DOUBLE_Q_BIT_SIZE - 1))) | (Q & -((d + (Q >> 1)) >> (DOUBLE_Q_BIT_SIZE - 1)));
  return d;
}

/*
 * Subtraction modulo q. Operands must be in the 0..q-1 range.
 */
S_DOUBLE_Q_SIZE sub(S_DOUBLE_Q_SIZE x, S_DOUBLE_Q_SIZE y)
{
  DOUBLE_Q_SIZE d;

  d = x - y;
  d += ((-Q) & ~-((d - ((Q >> 1) + 1)) >> (DOUBLE_Q_BIT_SIZE - 1))) | (Q & -((d + (Q >> 1)) >> (DOUBLE_Q_BIT_SIZE - 1)));
  return d;
}

/*
 * Division by 2 modulo q. Operand must be in the 0..q-1 range.
 */
S_DOUBLE_Q_SIZE rshift1(S_DOUBLE_Q_SIZE x)
{
  x += Q & -(x & 1);
  return (x >> 1);
}

/*
 * Montgomery multiplication modulo q. If we set R = 2^Q_BIT_SIZE mod q, then
 * this function computes: x * y / R mod q
 * Operands must be in the 0..q-1 range.
 */
S_DOUBLE_Q_SIZE montymul(S_DOUBLE_Q_SIZE x, S_DOUBLE_Q_SIZE y)
{
  DOUBLE_Q_SIZE x_, y_, z, w;

  x_ = (DOUBLE_Q_SIZE)addq(x);
  y_ = (DOUBLE_Q_SIZE)addq(y);

  z = x_ * y_;
  w = ((z * Q0I) & TWO_POWER_SIZE_Q_MINUS_ONE) * Q;

  z = (z + w) >> Q_BIT_SIZE;

  z -= Q;
  return reduce(z);
}

S_DOUBLE_Q_SIZE montymul_(S_DOUBLE_Q_SIZE x, S_DOUBLE_Q_SIZE y)
{
  DOUBLE_Q_SIZE z, w;

  z = x * y;
  w = ((z * Q0I) & TWO_POWER_SIZE_Q_MINUS_ONE) * Q;

  z = (z + w) >> Q_BIT_SIZE;

  z -= Q;
  z += Q & -(z >> (DOUBLE_Q_BIT_SIZE - 1));
  return z;
}

S_DOUBLE_Q_SIZE montymul2(S_DOUBLE_Q_SIZE x, S_DOUBLE_Q_SIZE y)
{
  return reduce((Q_SIZE)montymul_(montymul_(addq(x), addq(y)), R2));
}

/*************************************************
 * Name:        decompose
 *
 * Description: For finite field element a, compute high and low bits a0, a1 such
 *              that a mod^+ Q = a1*ALPHA + a0 with -ALPHA/2 < a0 <= ALPHA/2 except
 *              Assumes a to be standard representative.
 *
 * Arguments:   - S_Q_SIZE a: input element
 *              - S_Q_SIZE *a0: pointer to output element a0
 *
 * Returns a1.
 **************************************************/
S_Q_SIZE decompose(S_Q_SIZE *a0, S_Q_SIZE a)
{
  S_Q_SIZE a1, center, alpha = 2 * GAMMA2;

  // a -= Q * (a >> (Q_BIT_SIZE - 1));

  a1 = a >> (LOGGAMMA2 + 1);
  *a0 = a & (alpha - 1);

  center = ((alpha >> 1) - (*a0 + 1)) >> (Q_BIT_SIZE - 1);

  *a0 += alpha * center;
  a1 -= center;
  return a1;
}