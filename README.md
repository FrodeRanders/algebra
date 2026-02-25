# algebra

## Setting up a Python environment

```terminaloutput
➜  python -m venv .venv
➜  source .venv/bin/activate
(.venv) ➜  pip install maturin
Collecting maturin
  Using cached maturin-1.12.4-py3-none-macosx_10_12_x86_64.macosx_11_0_arm64.macosx_10_12_universal2.whl.metadata (16 kB)
Using cached maturin-1.12.4-py3-none-macosx_10_12_x86_64.macosx_11_0_arm64.macosx_10_12_universal2.whl (18.9 MB)
Installing collected packages: maturin
Successfully installed maturin-1.12.4

```

## Building the Algebra backend
```terminaloutput
(.venv) ➜  (cd algebrapy && maturin develop)
```

## Play with the Algebra package
```terminaloutput
(.venv) ➜  cd play
(.venv) ➜  python ghi.py
a = 3 (mod 7) b = 5 (mod 7)
a+b = 1 (mod 7)
a*b = 1 (mod 7)
inv(a) = 5 (mod 7)
a^6 = 1 (mod 7)
order(a) = 6
Fermat check OK

(.venv) ➜  cat rst.py 
import algebrapy as alg

K = alg.Fq(2, [1,1,0,1])      # 1 + x + x^3
a = K.elem([1,0,1])           # 1 + x^2
b = K.elem([0,1])             # x

print("a:", a)
print("b:", b)
print("a+b:", a + b)
print("a*b:", a * b)
print("b**7:", b**7)          # should be 1 in GF(2^3)^*
print("a/b:", a / b)
print("a**-1:", a**-1)
print("3 + a:", 3 + a)
print("2 * b:", 2 * b)        # note: in F2, 2 ≡ 0 => gives 0

(.venv) ➜  python rst.py 
a: 1 + x^2 in GF(2^3)
b: x in GF(2^3)
a+b: 1 + x + x^2 in GF(2^3)
a*b: 1 in GF(2^3)
b**7: 1 in GF(2^3)
a/b: 1 + x + x^2 in GF(2^3)
a**-1: x in GF(2^3)
3 + a: x^2 in GF(2^3)
2 * b: 0 in GF(2^3)

```

