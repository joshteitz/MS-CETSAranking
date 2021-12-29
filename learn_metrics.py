from metric_learn import ITML
from metric_learn import MMC
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

def learn_itml(tr, min_num = 10, use_neg_pairs = True):

  # if no training data
  if (tr is None):
    return None

  # if not enough training examples
  if ( len(tr["pairs_indices"]) < min_num ):
    return None
  
  # fit model to training data
  itml = ITML(preprocessor = tr["X"])
  if use_neg_pairs:
    itml.fit(tr["pairs_indices"], tr["y_pairs"])
  else:
    neg_index = tr["y_pairs"].index(-1)
    itml.fit(tr["pairs_indices"][:neg_index], tr["y_pairs"][:neg_index])

  return(itml.get_mahalanobis_matrix())
  
# def learn_mmc(tr, min_num = 10, use_neg_pairs = True):
#   
#   # if no training data
#   if (tr is None):
# 
#   # if not enough training examples
#   if ( len(tr["pairs_indices"]) < min_num ):
#     return None
#   
#   # fit model to training data  
#   mmc = MMC(preprocessor = tr["X"])
#   if use_neg_pairs:
#     mmc.fit(tr["pairs_indices"], tr["y_pairs"])
#   else:
#     neg_index = tr["y_pairs"].index(-1)
#     mmc.fit(tr["pairs_indices"][:neg_index], tr["y_pairs"][:neg_index])
#     
#   return(mmc.get_mahalanobis_matrix())
