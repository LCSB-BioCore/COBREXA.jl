"""
$(TYPEDEF)

This is the abstract type used to record the type of analysis performed.
"""
abstract type ResultType end

struct FluxBalanceAnalysis end

struct FluxVariabilityAnalysis end

struct MaxMinDrivingForceAnalysis end

struct Sampling end

