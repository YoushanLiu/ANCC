/*
This file is part of ANCC.
AFTAN is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
AFTAN is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
/*
 * Byte swap an array of length n of N byte integer*4 or real*4 elements.
 * A good compiler should unroll the inner loops. Letting the compiler do it
 * gives us portability.  Note that we might want to isolate the
 * cases N = 2, 4, 8 (and 16 for long double and perhaps long long)
 */
#define SWAP2(a, b) \
{               \
        (a) ^= (b); \
        (b) ^= (a); \
        (a) ^= (b); \
}

void swapn(unsigned char *b, int N, int n)
{

    int i, j;

    for (i = 0; i < (n) * (N); i += N)
        for (j = 0; j < N / 2; j++)
            SWAP2(b[i + j], b[i + N - j - 1]);

}
