import torch
from torch import nn
from copy import deepcopy
from tqdm import tqdm

class Original(torch.nn.Module):
    def __init__(self, n_steps, alpha, beta_bounds, tau_bounds, N0_bounds,
                ignore_diffs):
        super(Original, self).__init__()
        self.n_steps = n_steps
        self.alpha = alpha
        self.beta_bounds = beta_bounds
        self.tau_bounds = tau_bounds
        self.N0_bounds = N0_bounds
        self.ignore_diffs = ignore_diffs
        # declare and initialize trainable parameters
        self.beta = torch.nn.Parameter(torch.Tensor(n_steps))
        nn.init.uniform_(self.beta, self.beta_bounds[0], self.beta_bounds[1])
        self.tau = torch.nn.Parameter(torch.Tensor(1))
        nn.init.uniform_(self.tau, self.tau_bounds[0], self.tau_bounds[1])
        self.N0 = torch.nn.Parameter(torch.Tensor(1))
        nn.init.uniform_(self.N0, self.N0_bounds[0], self.N0_bounds[1])

    def clip_parameters(self):
        self.beta.data = torch.clamp(self.beta, self.beta_bounds[0], self.beta_bounds[1])
        self.tau.data = torch.clamp(self.tau, self.tau_bounds[0], self.tau_bounds[1])
        self.N0.data = torch.clamp(self.N0, self.N0_bounds[0], self.N0_bounds[1])

    def forward(self):
        self.clip_parameters()
        """
        expected_N = torch.zeros(self.n_steps, 1)
        expected_I = torch.zeros(self.n_steps, 1)
        for t in range(self.n_steps):
            if t == 0:
                expected_N[t] += self.N0 * (self.beta[t] + torch.exp(-self.tau))
                expected_I[t] += self.N0 * self.beta[t]
            else:
                expected_N[t] += expected_N[t - 1].clone() * (self.beta[t] + torch.exp(-self.tau))
                expected_I[t] += expected_N[t - 1].clone() * self.beta[t]
        return expected_I
        """
        expected_N = torch.zeros(self.n_steps)
        expected_I = torch.zeros(self.n_steps)
        for t in range(self.n_steps):
            expected_I[t] = self.beta[t] * self.N0 * torch.prod(self.beta[:t] + torch.exp(-self.tau))
            expected_N[t] = self.N0 * torch.prod(self.beta[:(t + 1)] + torch.exp(-self.tau))
        return expected_I

    def loss(self, expected_I, observed_I):
        beta_diffs = (self.beta[1:] - self.beta[:-1]) ** 2
        indices = torch.arange(len(beta_diffs))
        for d in self.ignore_diffs:
            beta_diffs = beta_diffs[indices != d]
        beta_smoother = self.alpha * torch.mean(beta_diffs)
        return torch.mean((expected_I - observed_I) ** 4) + beta_smoother

if __name__ == "__main__":
    # model parameters
    alpha = 2 ** 40
    beta_bounds = (0.1, 0.7)
    tau_bounds = (1/40, 1/4)
    N0_bounds = (0.01, 5)
    ignore_diffs = [26]
    observed_I = [1,0,3,3,2,2,1,8,10,5,7.31818181818182,12.5454545454545,16.7272727272727,18.39146567718,20.7823562152134,47.3375891568749,48.6599574589194,95.3254139490585,88.5025857099268,92.3009572961247,95.7255782506535,138.215173896726,187.339564908268,143.705768774486,116.271061378776,149.340765555692,136.22787207243,162.514923785303,176.519069949137,162.616184054277,199.470683669103,141.184155416023,163.06318449993,133.090672405605,254.012635658744,194.170223289247,194.391073504455,220.057711651108,191.101917670659,293.228206602399,245.045269893477,276.632712796882,281.307802551787,369.337022158363,338.928362500319,427.573695433027,465.274912414549,497.164933765501,445.874652320704]
    observed_I = torch.Tensor(observed_I)
    # optimization parameters
    lr = 10 ** -4
    max_rounds = 2 * (10 ** 5)
    scheduler_tolerance = 100
    stop_tolerance = 1000
    min_improve = 0.001

    model = Original(len(observed_I), alpha, beta_bounds, tau_bounds, N0_bounds, ignore_diffs)
    # instantiate optimization objects
    # optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    from radam import RAdam
    # from lookahead import Lookahead
    optimizer = Adam(model.parameters(), lr=lr) 
    # optimizer =  Lookahead(base_optim, k=5, alpha=0.5)
    # scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    #     optimizer, "min", patience=scheduler_tolerance, verbose=True
    # )
    # optimize loss
    losses = []
    last_improved = 0
    # train until convergence
    for r in tqdm(range(max_rounds)):
        torch.autograd.set_detect_anomaly(True)
        expected_I = model()
        loss = model.loss(expected_I, observed_I)
        losses.append(loss.item())
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        # scheduler.step(loss.item())
        # early stop if improvements stall
        if r > 0:
            improve = (min(losses[:-1]) - losses[-1]) / min(losses[:-1])
            if improve > min_improve:
                last_improved = 0
                best_model = deepcopy(model)
            elif last_improved >= stop_tolerance:
                break
            else:
                last_improved += 1
    import pdb; pdb.set_trace()