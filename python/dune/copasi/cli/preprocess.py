import os
import configparser
import sympy
from sympy.parsing.sympy_parser import parse_expr

def inject_jacobian(section,config):
  jac_section = section+'.jacobian'

  jac_ditc = {}
  for vari_str in config[section]:
    func_str = config[section][vari_str].replace('^','**')
    func_symb = parse_expr(func_str)

    for varj_str in config[section]:
      jac_key = 'd({})/d({})'.format(vari_str,varj_str)

      has_jac = False
      if config.has_section(jac_section):
        if jac_key in config[jac_section]:
          jac_str = str(config[jac_section][jac_key])
          has_jac = True
      
      if not has_jac:
        varj_symb = sympy.Symbol(varj_str)
        jac_symb = sympy.diff(func_symb,varj_symb)
        jac_str = str(jac_symb).replace('**','^')

      jac_ditc[jac_key] = jac_str
  config[jac_section] = jac_ditc
  return config

def preprocess_compartement(compartement,config):
  reaction_key = 'model.'+compartement+'.reaction'
  config = inject_jacobian(reaction_key,config)
  return config

def preprocess_config(args):

  if not os.path.isfile(args["config"]):
    raise IOError("Configuration file {} not found".format(args["config"]))

  config = configparser.ConfigParser()
  config.read(args['config'])
  

  for compartement in config['model.compartements']:
    preprocess_compartement(compartement, config)
  
  with open(args['new_config'], 'w') as configfile:
    config.write(configfile)