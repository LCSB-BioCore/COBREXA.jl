#! format: off
function __init__()
    py"""
    try:
       from equilibrator_api import ComponentContribution, Q_
       cc = ComponentContribution()
    
       def pygetdg0(fs, ph, ionic):
          cc.p_h = Q_(ph)  # set pH
          cc.ionic_strength = Q_(ionic)  # set I
          bals = []
          mags = []
          errs = []
          for f in fs:
             try:
                   rxn = cc.parse_reaction_formula(f)
                   isbal = rxn.is_balanced()
                   bals += [isbal]
                   mags += [cc.standard_dg(rxn).value.magnitude]
                   errs += [cc.standard_dg(rxn).error.magnitude]
             except:
                   bals += [False]
                   mags += [0.0]
                   errs += [0.0]
                   
          return bals, mags, errs
    
       def pygetdgprime(fs, ph, ionic):
          cc.p_h = Q_(ph)  # set pH
          cc.ionic_strength = Q_(ionic)  # set I
          bals = []
          mags = []
          errs = []
          for f in fs:
             try:
                   rxn = cc.parse_reaction_formula(f)
                   isbal = rxn.is_balanced()
                   bals += [isbal]
                   mags += [cc.standard_dg_prime(rxn).value.magnitude]
                   errs += [cc.standard_dg_prime(rxn).error.magnitude]
             except:
                   bals += [False]
                   mags += [0.0]
                   errs += [0.0]
                   
          return bals, mags, errs
          
       def pygetdgprimephys(fs, ph, ionic):
          cc.p_h = Q_(ph)  # set pH
          cc.ionic_strength = Q_(ionic)  # set I
          bals = []
          mags = []
          errs = []
          for f in fs:
             try:
                   rxn = cc.parse_reaction_formula(f)
                   isbal = rxn.is_balanced()
                   bals += [isbal]
                   mags += [cc.physiological_dg_prime(rxn).value.magnitude]
                   errs += [cc.physiological_dg_prime(rxn).error.magnitude]
             except:
                   bals += [False]
                   mags += [0.0]
                   errs += [0.0]
                   
          return bals, mags, errs
    
    except ImportError as e: # in case equilibrator is not installed.
       def pygetdg0(fs, ph, ionic):
          bals = [0.0]
          mags = [0.0]
          errs = [0.0]
          return bals, mags, errs
    
       def pygetdgprimephys(fs, ph, ionic):
          bals = [0.0]
          mags = [0.0]
          errs = [0.0]
          return bals, mags, errs
          
       def pygetdgprime(fs, ph, ionic):
          bals = [0.0]
          mags = [0.0]
          errs = [0.0]
          return bals, mags, errs
    """
end
