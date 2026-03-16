import sys, math
import numpy as np

sys.path.insert(0, '/Users/axel/axel/code/IKDH/tools')

from eval_jacobian import load_jacobian, evaluate_jacobian

rows = load_jacobian('/Users/axel/axel/code/IKDH/robots/gofa5.yaml')

q_deg = [0, 0, 0, 0, 0, 0]
q_rad = [math.radians(v) for v in q_deg]

J   = evaluate_jacobian(rows, q_rad)
det = np.linalg.det(J)

print(f"q = {q_deg} deg")
print(f"|det(J)| = {abs(det):.4e}")