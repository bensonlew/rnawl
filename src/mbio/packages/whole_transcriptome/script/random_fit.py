import logging

import matplotlib.pyplot as plt
import numpy as np
from keras import backend as K
from keras.applications import VGG16

logging.basicConfig(format='%(asctime)s\t%(name)s\t%(levelname)s : %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

model = VGG16(weights='imagenet', include_top=False)


size = 100
mu = 70.3
sigma = 5.25
rate = 0.1

arr = np.random.normal(mu, sigma, size)

def check(arr):
    if round(arr.mean(), 2) == mu and round(arr.std(), 2) == sigma:
        return True
    else:
        return False


while not check(arr):
    index = np.random.randint(0, size)






if __name__ == '__main__':
    pattern_dict = generate_patterns(['block1_conv1', 'block2_conv1', 'block3_conv1', 'block4_conv1'])
    for layer_name, results in pattern_dict.items():
        plt.figure(figsize=(24, 24))
        plt.title(layer_name)
        plt.imshow(results)
        plt.waitforbuttonpress()
        plt.close()
