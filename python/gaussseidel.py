"""
Descrição: código que soluciona qualquer equação do tipo Ax=b, para qualquer matriz A simétrica e positiva definida.
"""


def gaussseidel(A, b, erro):
    x0 = [0]*len(A)
    x = [None]*len(A)
    parada = False
    while not parada:
        for i in range(len(A)):
            x[i] = (b[i] - soma1(A, x, i) - soma2(A, x0, i, len(A)))/A[i][i]
        if norma(subtracao(b, produtoMatrizVetor(A, x0))) < erro:
            parada = True
        x0 = x
    print(f'A solução é: {x}')
    return x


def soma1(A, x, i):
    soma = 0
    for j in range(i):
        soma += A[i][j]*x[j]
    return soma


def soma2(A, x, i, n):
    soma = 0
    for j in range(i+1, n):
        soma += A[i][j]*x[j]
    return soma


def norma(v):

    """
    Descrição: calcula a norma euclidiana de qualquer vetor v;

    Entrada(s):
                i) v (list): vetor;
    
    Saída(s):
                i) modulo (float): norma euclidiana do vetor.
    """

    modulo = 0
    for i in v:
        modulo += i**2
    modulo = modulo**0.5
    return modulo


def subtracao(u, v):

    """
    Descrição: opera a subtração entre u e v. Por meio de zip, acessa as entradas de ambos os vetores
                e adiciona a diferença num outro vetor declarado, s;
    
    Entrada(s):
                i) u (list): vetor a ser subtraido (minuendo);
                ii) v (list): vetor subtraendo;
    
    Saída(s):
                i) s (list): vetor diferença.
    """

    s = list()
    for i, j in zip(u, v):
        s.append(i - j)
    return s


def produtoescalar(u, v):

    """
    Descrição: opera produto escalar entre u e v. Por meio de zip, acessa as entradas de ambos os vetores
                e adiciona o produto entre ambas entradas na variável pe;

    Entrada(s):
                i) u (list): vetor operador;
                ii) v (list): vetor operado;

    Saída(s):
                i) pe (float): produto escalar entre u e v.
    """

    pe = 0
    for i, j in zip(u, v):
        pe += i*j
    return pe


def produtoMatrizVetor(A, v):

    """
    Descrição: opera o produto entre uma matriz e um vetor, nessa ordem;

    Entrada(s):
                i) A (list): matriz;
                ii) v (list): vetor;
    
    Saída(s):
                i) p (list): vetor resultado da operação.
    """

    p = list()
    for i in A:
        p.append(produtoescalar(i, v))
    return p

A = [[4, -1, 1], [-1, 5, 2], [1, 2, 9]]  # Ser SPD!
b = [2, 10, -1]
gaussseidel(A, b, 1e-6)
