import os

for seed in range(100):

    vis_command = "java StarTravellerVis -exec \"/home/dawid/TC/StarTraveller/./StarTravellerExtended\" -seed "
    vis_command = vis_command + str(seed)
    os.system(vis_command)
    print("\n")