#include "../inc/common.h"

namespace GridTool::COMMON
{
    Scalar relaxation(Scalar a, Scalar b, Scalar x)
    {
        return (1.0 - x) * a + x * b;
    }
}
