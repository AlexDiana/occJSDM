# Called by sjsdm-gradient-noise.R with the original native loss and fixed weights.
import csv
import time
import numpy as np
import torch

m = diagnostic_model
x = torch.as_tensor(np.asarray(diagnostic_x), dtype=torch.float64)
y = torch.as_tensor(np.asarray(diagnostic_y), dtype=torch.float64)
parameters = list(m.env.parameters()) + [m.sigma]
direction = torch.as_tensor(np.asarray(diagnostic_direction).reshape(-1), dtype=torch.float64)
assert m._loss_function.__name__ == "torch_tmp"
assert m.df == 2 and str(m.device) == "cpu"
assert len(parameters) == 2 and parameters[0].shape == (10, 3)
assert parameters[1].shape == (10, 2)
assert all(t.dtype == torch.float64 for t in parameters)
torch.set_num_threads(1)
started = time.monotonic()
rows = []
for sampling, repeats in [(2000, 100), (20000, 40)]:
    for repeat in range(repeats):
        seed = int(diagnostic_seed) + sampling + repeat
        torch.manual_seed(seed)
        loss = m._loss_function(m.env(x), y, m.sigma, 100, sampling, 2, 1.0, "cpu", torch.float64).sum()
        penalised_loss = loss + 100 * 0.0001 / 2 * sum((p*p).sum() for p in parameters)
        gradients = torch.autograd.grad(penalised_loss, parameters)
        gradient = torch.cat([g.reshape(-1) for g in gradients])
        row = [sampling, repeat + 1, seed, float(penalised_loss.detach()),
               float(torch.dot(gradient, direction))] + gradient.detach().tolist()
        rows.append(row)
    print("Native gradient experiment", sampling, "draws finished in", time.monotonic()-started, "seconds", flush=True)
with open(diagnostic_output, "w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["sampling", "repeat", "seed", "penalised_loss", "directional_gradient"] +
                    ["gradient_" + str(k + 1) for k in range(50)])
    writer.writerows(rows)
