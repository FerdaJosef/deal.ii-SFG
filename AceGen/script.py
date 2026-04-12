import re
import sys

def transform_acegen_file(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # This regex is much looser:
    # 'void\s+equation' -> matches 'void equation'
    # '\(.*?\)' -> matches EVERYTHING inside the first set of parentheses, including newlines
    # '\s*\{' -> matches any whitespace and the opening bracket
    pattern = r"void\s+equation\s*\(.*?\)\s*\{"

    new_header = """inline void equation_wrapper(
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

    # Check if we even find a match before trying to write
    if re.search(pattern, content, flags=re.DOTALL):
        new_content = re.sub(pattern, new_header, content, count=1, flags=re.DOTALL)
        with open(filename, 'w') as f:
            f.write(new_content)
        print(f"Successfully transformed {filename}")
    else:
        print("Error: Could not find the function 'void equation' in the file.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        transform_acegen_file(sys.argv[1])
    else:
        print("Usage: python3 script.py equation0.h")