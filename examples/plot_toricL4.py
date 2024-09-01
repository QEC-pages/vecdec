import matplotlib.pyplot as plt
import numpy as np

# Data
physical_error_rate = [0.005, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1,0.108,0.12, 0.2, 0.3, 0.4]
logical_error_rate_ML = [0.00076, 0.00383, 0.0102, 0.0186, 0.03, 0.117, 0.3704,0.417333,0.469 ,0.678, 0.6925, 0.704375]
logical_error_rate_mode0 = [0.00056, 0.00237, 0.00612, 0.01108, 0.01794, 0.07603, 0.27722,0.31507 ,0.36712 ,0.62878, 0.73311, 0.75083]


plt.figure(figsize=(10, 6))

# double-log
plt.plot(physical_error_rate, logical_error_rate_ML, '-^', markersize=8, linewidth=2, label='ML')
plt.plot(physical_error_rate, logical_error_rate_mode0, '-o', markersize=8, linewidth=2, label='Mode 0')

#
plt.title('Logical Error Rate vs Physical Error Rate', fontsize=16)
plt.xlabel('Physical Error Rate', fontsize=14)
plt.ylabel('Logical Error Rate', fontsize=14)

# grid on
plt.grid(True, which="both", ls="-", alpha=0.2)

# 
plt.xlim(0.004, 0.5)
plt.ylim(0.0001, 1)

#
plt.legend(loc='best')

# 
plt.savefig('error_rate_plot.png', dpi=300, bbox_inches='tight')

# 
plt.show()