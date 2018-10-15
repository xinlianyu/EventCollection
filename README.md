# EventCollection
generate event collection in C++


We study optimal routing policy problems in stochastic time-dependent (STD) networks, where link travel times are modeled as random variables with time-dependent distributions. A routing policy is defined as a decision rule that specifies what node to take next at each
decision node based on realized link travel times and the current time.

The concept of information set is introduced to represent the traveler's knowledge about the network. An information set is composed of support points that are consistent with the link travel times revealed so far. A "support point" is defined as a distinct value that a discrete random variable can take or a distinct vector of values that a discrete random vector can take, depending on the context. Thus a probability mass function (PMF) of a random variable (vector) is a combination of support points and the associated probabilities. When the information set becomes a singleton, the network becomes deterministic.

Define  the event collection EV as the set of support point candidates after collecting information at time t. As we collect more information (i.e. t increases), the size of EV remains the same or decreases.
