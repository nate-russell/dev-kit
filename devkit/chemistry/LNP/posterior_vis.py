import torch
import gpytorch
from botorch.models import SingleTaskGP
from botorch.test_functions import Branin
from botorch.utils import draw_sobol_samples

# Define a test function (Branin function in this case)
test_function = Branin()

# Generate training data using Sobol samples
train_X = draw_sobol_samples(test_function.bounds,n=5,q=1)
train_Y = test_function(train_X)

# Define the Gaussian Process model using GPyTorch
class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_X, train_Y):
        likelihood = gpytorch.likelihoods.GaussianLikelihood()
        super(ExactGPModel, self).__init__(train_X, train_Y, likelihood)

    def forward(self, x):
        mean = self.likelihood(self.likelihood(x))
        return gpytorch.distributions.MultivariateNormal(mean, self.covariance_matrix)

# Create the GPyTorch model
gpytorch_model = ExactGPModel(train_X, train_Y)

# Convert GPyTorch model to Botorch model
botorch_model = SingleTaskGP(gpytorch_model)

# Visualize the posterior distribution
with torch.no_grad():
    posterior_samples = gpytorch_model.posterior(train_X).sample(sample_shape=torch.Size([10]))

# Plot the training data and posterior samples
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
plt.scatter(train_X.numpy(), train_Y.numpy(), color='red', marker='x', label='Training Data')

for i in range(posterior_samples.shape[0]):
    plt.plot(train_X.numpy(), posterior_samples[i].numpy(), linestyle='-', alpha=0.3, label='Posterior Sample')

plt.title('Posterior Distribution of GP')
plt.xlabel('Input')
plt.ylabel('Output')
plt.legend()
plt.show()
