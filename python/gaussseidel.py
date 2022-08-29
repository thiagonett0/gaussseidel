"""
Descrição: código que soluciona qualquer equação do tipo Ax=b, para qualquer matriz A simétrica e positiva definida.
"""


import sys
sys.path.append( 'path' )

import linalg as la


def soma1(A, x, i):

    """
    Descrição: calcula o primeiro somatório apresentado na última expressão do
                arquivo readme.md;
    
    Entrada(s):
                i) A (list): matriz do sistema;
                ii) x (list): vetor iterando, no passo posterior;
                iii) i (int): i-ésima entrada do vetor x;
    
    Saída(s):
                i) soma (float): resultado do somatório.
    """

    soma = 0
    for j in range(i):
        soma += A[i][j]*x[j]
    return soma


def soma2(A, x, i):

    """
    Descrição: calcula o segundo somatório apresentado na última expressão do
                arquivo readme.md;
    
    Entrada(s):
                i) A (list): matriz do sistema;
                ii) x (list): vetor iterando, no passo anterior;
                iii) i (int): i-ésima entrada do vetor x;
    
    Saída(s):
                i) soma (float): resultado do somatório.
    """

    soma = 0
    for j in range(i+1, len(A)):
        soma += A[i][j]*x[j]
    return soma


def gaussseidel(A, b, erro):

    """
    Descrição: calcula a solução aproximada x de um sistema Ax=b, para A simétrica e positiva definida. 
            O algoritmo segue conforme o discutido em Teorema 10.2.1, página 512 de Matrix Computations, 
            Golub e Van Loan, com algumas observações:

                i) arbitra-se como estimativa inicial o vetor nulo;
                ii) adota-se um critério de parada conservador, como segue no código;
                iii) o funcionamento do código é análogo a qualquer processo iterativo, i.e.:

                        iii.a) após estabelecidos a estimativa inicial, a tolerância e o critério de parada, calcula
                                a solução no próximo passo por meio da expressão pré determinada;
                        iii.b) verifica-se se a solução satisfaz o critério de parada. Se sim, altera a variável 
                                booleana parada;
                        iii.c) a variável que guardava o passo anterior é alterada de forma a guarda a solução no
                                no passo atual;
                        iii.d) repete até a variável parada ser alterada;

    Entrada(s):
                i) A (list): matriz simétrica e positiva definida;
                ii) b (list): vetor independente do sistema;
                iii) erro (float): tolerância de erro aceitável;
    
    Saída(s):
                i) x (list): solução aproximada do sistema.
    """

    x0 = [0]*len(A)
    x = [None]*len(A)
    parada = False
    while not parada:
        for i in range(len(A)):
            x[i] = (b[i] - soma1(A, x, i) - soma2(A, x0, i))/A[i][i]
        if la.euclidianNorm(la.vectorSum(b, la.vectorScalar(-1, la.matrixVector(A, x0)))) < erro:
            parada = True
        x0 = x
    print(f'A solução é: {x}')
    return x