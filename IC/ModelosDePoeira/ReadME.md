Poeira:

- Separada em componentes de emissão e absorção, sendo a emissão separada em duas componentes, uma de poeira quente e outra de poeira fria e difusa (de nuvens moleculares, é da componente nebular, gerada por um código de fotoionização Cloudy).
- Transforma luz do azul/ultravioleta em infravermelho próximo
- Modelos usam poeira com silicatos, grafite e PAH (hidrocarbonetos aromáticos policíclicos), um tipo de poeira orgânica que reemite e a luz por radiação térmica que depende do tamanho da partícula. Maiores -> menor comprimento de onda, mais quente.
- A emissão de PAHs pequenos geralmente acontece por fóton único, então é independente da intensidade da luz estelar (U, relacionado com o parâmetro Umin do bagpipes), principalmente para U < 1000.
- Modelos atuais de poeira tendem a subestimar a emissão na faixa dos 4.5 micrômetros
- QPahs são relevantes apenas para lambdas altos, não tendo efeitos no visível, assim como umin e gamma, que são variáveis que definem a intensidade da radiação que aquece essa poeira

- Transmitância: log10(T0(aᵢ, λ)) = −0.4 * 𝜖 * A_V * k(λ) / R_V (Estrelas jovens têm o parâmetro 𝜖 (ou eta), uma constante que diz que a poeira ali é maior do que o normal, por ter passado por uma nuvem molecular, portanto populações estelares com: aᵢ (idade da estrela) < a_BC (tempo que a nuvem molecular fica no entorno) têm o parâmetro 𝜖 na equação.)
- k(λ) é a curva de atenuação e R_V é a razão da atenuação total pelo avermelhamento. Estes valores são definidos por cada modelo, sendo extrapolada para comprimentos de onda menores usando um ajuste exponencial pra curva de poeira no ultravioleta próximo em Calzetti e Cardelli.

- Calzetti: Atenua o UV e o óptico menos que, por exemplo, a poeira da Via Láctea (Cardelli). A atenuação em direção às regiões HII é 2-3 vezes maior do que para as estrelas mais velhas. Foi desenvolvida para galáxias locais com formação estelar intensa e no Bagpipes é usada em modelos para galáxias com formação estelar ativa.

- Cardelli: Atenua mais fortemente no UV, principalmente em 2175 Å, pois foi desenvolvido para modelar a poeira da Via Láctea, que tem uma propriedade relacionada à poeira com grãos de grafite. É útil no bagpipes para estudar galáxias espirais semelhantes à Via Láctea.

- Charlot & Fall (CF00): É mais flexível, com k(λ)/R_V = (λ/5500Å)^(-n), onde podemos variar o valor de "n" para obter curvas de atenuação bem diferentes, possibilitando usar modelos de poeira raros ou até desconhecidos.

- Salim: Segue a equação (k_lambda,mod = k_lambda,Cal * (lambda/5500Å)^delta + B * D_lambda * R_V,mod / R_V,Cal) que adiciona uma "corcova UV" mais proeminente quanto maior for o valor de B e mais inclinada (inclinação da lei exponencial atenuando mais o azul e menos o vermelho, por exemplo) conforme o valor de delta. Essa "corcova" causa uma queda acentuada na luminosidade observada na região do ultravioleta próximo.
Como as curvas de atenuação podem ser mais ou menos íngremes com o modelo de Salim, ele pode reproduzir o modelo de poeira de galáxias de alto redshift, alta opacidade óptica, baixa metalicidade, quiescentes e com formação ativa de estrelas (como a da SMC).


Síntese de População Estelar:

- Supõe que a luz da galáxia é formada por um grupo de estrelas com a mesma idade e metalicidade de um lugar determinado
- Desta suposição, forma a isócrona destas estrelas e pega o espectro das estrelas usando bibliotecas com abordagens teóricas e outras com abordagens com dados reais
- Dá uma maior importância para TPAGBs (Thermally Pulsing Asymptotic Giant Branch) - estrelas gigantes vermelhas no fim da vida com massa intermediária - são muito luminosas no infravermelho e ejetam poeira e elementos químicos enriquecidos, contribuindo muito para a modelagem da luz de galáxias, principalmente se forem formadas recentemente
- Age-metallicity Degenerescence: população estelar velha e pobre em metais pode ter uma cor integrada muito parecida com a de uma população jovem e rica em metais. Para contornar este problema pode ser medida a força de linhas de absorção específicas -> índices espectrais; cobertura maior de comprimentos de onda

Multinest (Nested Sampling):

- Vasculha por um mapa de nível com "montanhas" e "vales", calculando pontos com maior verossimilhança e formando um caminho até o topo da maior montanha que puder achar dentro dos limites dos prioris definidos, onde cada ponto define um modelo diferente do bagpipes
- Usa outro tipo de algoritmo que guarda os pontos de menor verossimilhança que normalmente seriam descartados para poder ter um mapa final muito melhor (melhor resolução da densidade de pontos), calculando com maior precisão a evidência baysiana (nota final do modelo).
- Quanto maior o número de parâmetros para ajustar, maior a dimensionalidade do problema e maior o custo computacional, então para o ajuste de espectros é esperado que alguns parâmetros específicos possam ser definidos de acordo com a literatura ao invés de serem ajustados.
