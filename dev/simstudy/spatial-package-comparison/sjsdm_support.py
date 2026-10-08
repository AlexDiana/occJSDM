"""Study-only checks and serialisation around the pinned native sjSDM model."""
import numpy as np
import torch

def parameters(model):
    return list(model.env.parameters()) + ([] if model.spatial is None else list(model.spatial.parameters()))

def check_diagonal(model):
    assert model.sigma.shape[0] == model.df
    assert torch.equal(model.sigma, torch.eye(model.df, dtype=torch.float64))
    assert not model.sigma.requires_grad
    expected = {id(p) for p in parameters(model)}
    actual = {id(p) for g in model.optimizer.param_groups for p in g['params']}
    assert expected == actual and id(model.sigma) not in actual
    assert all(p.dtype == torch.float64 and str(p.device) == 'cpu' for p in parameters(model))
    return True

def restart_optimizer(model, lr, weight_decay):
    check_diagonal(model)
    model.optimizer = torch.optim.RMSprop(parameters(model), lr=lr, alpha=.99,
        eps=1e-8, momentum=.1, centered=False, weight_decay=weight_decay)
    model.useSched = False
    check_diagonal(model)

def network_spec(network):
    if network is None:
        return []
    layers = []
    for layer in network:
        if isinstance(layer, torch.nn.Linear):
            layers.append(dict(kind='linear', weight=layer.weight.detach().numpy().copy(),
                bias=None if layer.bias is None else layer.bias.detach().numpy().copy()))
        elif isinstance(layer, torch.nn.SELU):
            layers.append(dict(kind='selu'))
        else:
            raise ValueError('Uncertified network layer: '+str(type(layer)))
    return layers

def snapshot(model):
    check_diagonal(model)
    return dict(environment=network_spec(model.env), spatial=network_spec(model.spatial),
                sigma=model.sigma.detach().numpy().copy(), df=model.df)

def exact_objective(model, X, Y, SP, l1, l2, weight_decay, nodes=60):
    check_diagonal(model)
    X = torch.as_tensor(X,dtype=torch.float64); Y=torch.as_tensor(Y,dtype=torch.float64)
    mu = model.env(X)
    if model.spatial is not None:
        mu = mu + model.spatial(torch.as_tensor(SP,dtype=torch.float64))
    z,w=np.polynomial.hermite_e.hermegauss(int(nodes)); w=w/np.sqrt(2*np.pi)
    z=torch.as_tensor(z,dtype=torch.float64); w=torch.as_tensor(w,dtype=torch.float64)
    p=(torch.sigmoid(1.70169*(mu[...,None]+z))*.999999+.0000005) @ w
    nll=-(Y*torch.log(p)+(1-Y)*torch.log1p(-p)).sum()/len(X)
    params=parameters(model)
    spatial=[] if model.spatial is None else list(model.spatial.parameters())
    penalized_ids={id(p) for p in spatial if p.ndim>1}
    # Native spatial elastic net penalises matrices, not biases. RMSprop weight
    # decay acts on every optimised parameter, including biases and intercepts.
    elastic=sum((l1*p.abs().sum()+l2*(p*p).sum() for p in spatial if p.ndim>1),torch.tensor(0.,dtype=torch.float64))
    native_penalty=sum((f() for f in model.losses),torch.tensor(0.,dtype=torch.float64))
    assert abs(float(native_penalty.detach()-elastic.detach()))<1e-9
    smooth=nll+sum((l2*(p*p).sum() for p in spatial if p.ndim>1),torch.tensor(0.,dtype=torch.float64))
    decay=sum((weight_decay*.5*(p*p).sum() for p in params),torch.tensor(0.,dtype=torch.float64))
    grads=torch.autograd.grad(smooth+decay,params)
    residuals=[]; spatial_gradient=[]
    for param,g in zip(params,grads):
        if id(param) in penalized_ids:
            residual=torch.where(param!=0,g+l1*torch.sign(param),torch.sign(g)*torch.clamp(g.abs()-l1,min=0))
            spatial_gradient.append(float(residual.abs().max().detach()))
        else:
            residual=g
        residuals.append(float(residual.abs().max().detach()))
    return dict(objective=float((nll+elastic+decay).detach()),nll=float(nll.detach()),
        elastic_penalty=float(elastic.detach()),weight_decay_penalty=float(decay.detach()),
        max_gradient=max(residuals),spatial_gradient=max(spatial_gradient,default=0.),
        probabilities=p.detach().numpy(),quadrature_nodes=int(nodes),
        normalisation='negative summed Bernoulli log likelihood / training sites + native elastic net + 0.5 weight_decay sum(all parameters squared)',
        stationarity='ordinary smooth gradient; L1 signed derivative for nonzero entries and distance to subgradient interval at exact zeros')

def certify_training_marginal(model, seed=9871):
    torch.manual_seed(int(seed))
    mu=torch.tensor([[-4.],[-1.],[0.],[1.],[4.]],dtype=torch.float64)
    batches=[]
    for _ in range(20):
        loss=model._loss_function(mu,torch.ones_like(mu),torch.eye(1,dtype=torch.float64),5,20000,1,1.70169,'cpu',torch.float64)
        batches.append(torch.exp(-loss).detach().numpy())
    a=np.array(batches)
    return dict(mu=mu[:,0].numpy(),mean=a.mean(axis=0),se=a.std(axis=0,ddof=1)/np.sqrt(len(a)),
                batches=20,draws_per_batch=20000,seed=int(seed))

def rng_state():
    import random
    return dict(torch=torch.get_rng_state().numpy(),numpy=np.random.get_state(),python=random.getstate())
