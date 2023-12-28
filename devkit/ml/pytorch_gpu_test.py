"""
figure out what GPU(s) are available and their memory
"""
import torch


def bytes_to_GB(b):
    gb = b / (1024 ** 3)
    return f"{gb:0.02f} GB"

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Current device: {device}')
n_devices = torch.cuda.device_count()
print(f"# of CUDA Devices = {n_devices}")
curent_device = torch.cuda.current_device()
print(f"Current CUDA Device = {curent_device}")

spacer = '-'*75
for i in range(n_devices):
    print(spacer)
    name = torch.cuda.get_device_name(i)
    total_memory = torch.cuda.get_device_properties(i).total_memory
    print(f"Device {i} Name: {name} Total Memory: {bytes_to_GB(total_memory)}")
    print(f'Allocated: {bytes_to_GB(torch.cuda.memory_allocated(0))}')
    print(f'Reserved: {bytes_to_GB(torch.cuda.memory_reserved(0))}')



