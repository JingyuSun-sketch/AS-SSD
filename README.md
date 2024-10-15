# AS-SSD
## Abstract
In high-speed network environments, super-spreaders are defined as hosts or devices with a large number of connections. Accurate super spreader detection is crucial for many applications such as network monitoring, security analysis, and traffic management. Invertible algorithms based on sketches have received extensive attention due to their excellent memory efficiency and the ability to directly recover the super spreader identity from the internal structure. According to application interests, the packets sent or received from the same host or device are abstracted as a flow. The flow distribution in real-world networks is highly skewed. Most flows are small, and only a few flows are large. However, existing solutions can hardly adapt to the highly skewed flow distribution efficiently, resulting in low memory efficiency. Therefore, this paper designs an adaptive sampling based super spreader detection algorithm AS-SSD. The algorithm proposes an adaptive sampling strategy based on register sharing, which makes up for the above shortcomings. AS-SSD first maps the arriving flow elements to a register array. Small flows only consume a few registers, while larger flows consume more, which adapts to the skewed flow distribution. Then, AS-SSD uses the adaptive sampling strategy to dynamically adjust the probability for sampling the flow elements and reduces the register updates caused by large flows. Therefore, it can avoid excessive memory occupation of large flows, thereby further improving the utilization efficiency of memory resources under the premise of ensuring accuracy. Experimental evaluation shows that compared with the previous work, AS-SSD shows higher detection accuracy in the super spreader detection task while maintaining high throughput processing capability. Compared with the most advanced algorithms, it can increase the F1 value by more than 0.609.

## About this repo
For AS-SSD, FreeRS, SpreadSketch, ExtendedSketch+, FlowFight, the compiling command is as follows (“***“ represents the name of the main algorithm code file):
```shell
g++ ***.cpp MurmurHash3.cpp -m64 -O3
```

For ExtendedSketch, the compiling command is as follows:
```shell
g++ ExtendedSketch.cpp MurmurHash3.cpp mini-gmp.c -m64 -O3
```



On win10 with MinGW, the compiling commands may need to add " -static-libstdc++ -static-libgcc".
