
def default_parameters():
    return {
        "enable_profiling": False,
        "enable_factorization": False,  # True, # Fails for hyperelasticity demo in dolfin, needs debugging
        "max_registers": 1024,  # 8 B * 1024 = 8 KB # TODO: Tune this for something very complex
        "score_threshold": 3,  # TODO: Scoring is work in progress and this will change meaning later
    }
