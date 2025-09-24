# Qlapoti

Code accompaining the paper "Qlapoti: simple and efficient translation of
quaternion ideal to isogenies" by Giacomo Borin, Maria Corte-Real Santos, Jonathan Komada Eriksen, Riccardo Invernizzi, Marzio Mula, Sina Schaeffler, and Frederik Vercauteren.

All code in this folder in the repository is writen in `Python` and [`SageMath`](https://www.sagemath.org/), and is made publicly available under the MIT license. The dependencies are:
- Python 3.11.1
- SageMath version 10.7

The main Qlapoti algorithm is implemented in `norm_eq.py`.

## Ideal to isogeny

The ideal to isogeny algorithm is implemented in
`ideal_to_isogeny/ideal_to_isogeny_qlapoti.py`. To verify the timings obtained
in the paper move to the folder `ideal_to_isogeny` and run `sage --python -O
compare.py`. Different security levels can be chosen inside the file.

## PRISM

The ideal to isogeny algorithm is integrated with the [sage implementation of
PRISM](https://github.com/KULeuven-COSIC/PRISM/tree/main/sage_implementation), in the folder `applications/PRISM`. To verify the timings obtained in
the paper go to that folder and run `sage --python -O test_prism.py`. Different
security levels can be chosen inside the file.

## SQISign

The ideal to isogeny algorithm is integrated with an unofficial proof of
concept implementation of SQIsign, which unfortunately is not provied as
the authors did not want it to be publically shared. The timings were 
measure in a similar manner as PRISM, except that the code currently only 
supports NIST level I.
