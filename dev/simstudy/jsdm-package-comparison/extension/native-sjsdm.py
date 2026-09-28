"""Expose the covariance already inverted by pinned sjSDM's native se routine.

The scoped observer does not alter torch.inverse's arguments or result. It is
restored even after an error. No package source or model weights are changed.
Only use in the serial study worker, with the frozen PyTorch CPU backend.
"""
from unittest.mock import patch
import numpy as np
import torch


def native_covariance(model, x, y, sampling, seed):
    matrices = []
    original = torch.inverse

    def observe(matrix, *args, **kwargs):
        result = original(matrix, *args, **kwargs)
        matrices.append((matrix.detach().cpu().numpy().copy(),
                         result.detach().cpu().numpy().copy()))
        return result

    before_beta = model.env_weights[0].copy()
    before_loading = model.get_sigma.copy()
    torch.manual_seed(int(seed))
    with patch.object(torch, "inverse", observe):
        se = np.asarray(model.se(np.asarray(x), np.asarray(y),
            batch_size=len(x), parallel=0, sampling=int(sampling), verbose=False))
    assert torch.inverse is original
    assert np.array_equal(before_beta, model.env_weights[0])
    assert np.array_equal(before_loading, model.get_sigma)
    assert len(matrices) == y.shape[1]
    covariance = np.asarray([item[1] for item in matrices])
    np.testing.assert_allclose(se**2, np.diagonal(covariance, axis1=1, axis2=2),
                               rtol=1e-12, atol=1e-12)
    return dict(se=se, covariance=covariance,
                regularized_hessian=np.asarray([item[0] for item in matrices]))
