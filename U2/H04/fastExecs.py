import numpy as np

# Caminho para o arquivo NPZ
path = r'C:\Users\pedro\Documents\UFRN\2024.2\Comunicações Móveis\Projetos\HandsOnCM\U2\H04\TBS_data.npz'

# Carregar o arquivo NPZ
tbs_data = np.load(path)

# Verificar as chaves disponíveis no arquivo NPZ
print(tbs_data.files)

# Acessar a matriz 'TBS'
tbs_matrix = tbs_data['TBS']

# Imprimir o elemento na posição [1][1]
print(len(tbs_matrix[0]))