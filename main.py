import os
import time

total = 0.0
N = 21 # number of cases (seeds) to test

for seed in range(1,N):

    vis_command = "java StarTravellerVis -exec \"/home/dawid/TC/StarTraveller/./StarTravellerCombinedMethods\" -seed "
    vis_command = vis_command + str(seed)

    start_time = time.time()
    output = os.popen(vis_command).readlines()
    finish_time = time.time()
    time_elapsed = finish_time - start_time

    print("Case " + str(seed) + " time: " + str(time_elapsed) + " score: " + output[-1], end="")

    total = total + float(output[-1])

mean = total/N
print("Mean score: " + str(mean))
