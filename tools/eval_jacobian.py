import math
import numpy as np
import yaml


def load_jacobian(yaml_path):
    with open(yaml_path) as f:
        return yaml.safe_load(f)['jacobian']['rows']


def evaluate_jacobian(rows, q):
    """Return the 6x6 geometric Jacobian (numpy array) at joint angles q (radians)."""
    env = {f's{i+1}': math.sin(q[i]) for i in range(6)}
    env.update({f'c{i+1}': math.cos(q[i]) for i in range(6)})
    env.update({f'q{i+1}': q[i] for i in range(6)})
    env['sin'] = math.sin
    env['cos'] = math.cos
    return np.array([[eval(e, {"__builtins__": {}}, env) for e in row] for row in rows])


if __name__ == '__main__':
    import sys
    yaml_path = sys.argv[1] if len(sys.argv) > 1 else 'robots/gofa5.yaml'
    rows = load_jacobian(yaml_path)
    q = [math.radians(float(v)) for v in sys.argv[2:]] if len(sys.argv) > 2 else [0.0] * 6
    J = evaluate_jacobian(rows, q)
    print(np.array2string(J, precision=4, suppress_small=True, separator=', '))
    print(f"det={np.linalg.det(J):.4e}  manip={math.sqrt(max(0, np.linalg.det(J @ J.T))):.4e}")