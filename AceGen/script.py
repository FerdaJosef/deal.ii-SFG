import re
import sys
import os

def transform_acegen_file(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # --------------------------------------------------
    # 1. Remove sms.h include
    # --------------------------------------------------
    content = re.sub(r'#include\s+"sms\.h".*\n', '', content)

    # --------------------------------------------------
    # 2. Replace Power(x,y) → std::pow(x,y)
    # --------------------------------------------------
    content = re.sub(r'\bPower\s*\(', 'std::pow(', content)

    # --------------------------------------------------
    # 3. Replace function header (ANY function name)
    # --------------------------------------------------
    pattern = r"void\s+\w+\s*\(.*?\)\s*\{"

    new_header = """template <int dim, int n>
inline void equation(
    std::vector<double> &v,
    const Vector<double> &U,
    const Vector<double> &U0,
    const std::vector<Tensor<1,dim>> &GradU,
    Vector<double> &dPsiDu,
    std::vector<Tensor<1,dim>> &dPsidGradU,
    FullMatrix<double> &dPsiDu2,
    std::vector<std::vector<Tensor<1,dim>>> &dPsidUdGradU,
    std::vector<std::vector<Tensor<2,dim>>> &dPsidGradU2,
    double (*dt))
{"""

    if re.search(pattern, content, flags=re.DOTALL):
        content = re.sub(pattern, new_header, content, count=1, flags=re.DOTALL)
    else:
        print("Warning: function header not found!")

    # --------------------------------------------------
    # 4. Fix possible AceGen indexing bug: v[5,3] → v[5][3]
    # --------------------------------------------------
    content = re.sub(r'(\w+)\[(\d+),\s*(\d+)\]', r'\1[\2][\3]', content)

    # --------------------------------------------------
    # 5. Build header structure
    # --------------------------------------------------
    header = """#pragma once

#include <vector>
#include <cmath>

#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

using namespace dealii;

"""

    content = header + content

    # --------------------------------------------------
    # 6. Write to NEW .h file
    # --------------------------------------------------
    base, _ = os.path.splitext(filename)
    new_filename = base + ".h"

    with open(new_filename, 'w') as f:
        f.write(content)

    print(f"Processed {filename} → {new_filename}")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        transform_acegen_file(sys.argv[1])
    else:
        print("Usage: python3 script.py equation0.c")