# Redesenho do FamFlow: resumo para revisão

_Leitura de uns 15 minutos. As seções 1 a 4 explicam a proposta. A seção 5 mostra, de forma geral, como o código será reorganizado. **A seção 6 traz as decisões que dependem de você.**_

---

## 1. A ideia em um parágrafo

O FamFlow deve ajudar o pesquisador a descobrir **como uma família de proteínas está organizada de verdade**, e não só a parte parecida com a sequência de partida. A partir de um conjunto de homólogos, ele encontra **todas as comunidades** da rede de similaridade de sequências (SSN). Para cada comunidade, informa o **tamanho real** e as características: taxonomia, arquitetura de domínios e vizinhança genômica. **Sequências órfãs**, as mais divergentes, são tratadas como candidatas, não como ruído. Depois disso, o FamFlow permite **testar** se uma comunidade, ou qualquer grupo que você definir, é de fato uma **família**. Para isso, ele constrói um HMM de perfil para cada grupo e mede se esses HMMs **separam** os grupos de origem.

O julgamento biológico continua sendo do pesquisador: o que entra, quais cortes usar, quais marcadores contam. O papel do FamFlow é deixar cada uma dessas escolhas **explícita, registrada e reversível**, e permitir parar em qualquer etapa com um resultado utilizável.

---

## 2. O que a ferramenta faz

O fluxo continua o mesmo de hoje. Cada etapa lê o resultado da anterior, grava o próprio resultado e um relatório curto, e pode ser o ponto final da análise.

```
Semente (uma sequência, um conjunto, um HMM pronto, ou uma sequência fatiada em domínios)
  │
  ├─ search      busca por HMM no banco de referência
  ├─ inspect     gráficos de score × cobertura do perfil, antes de qualquer corte
  ├─ collect     guarda os hits que passam nos filtros (só o envelope do domínio ou a proteína inteira)
  ├─ cluster     remove sequências quase idênticas, mas lembra quem cada representante representa
  ├─ matrix      similaridade todos contra todos entre os representantes
  ├─ ssn         rede + comunidades + órfãs                         ← bom ponto de parada
  ├─ annotate    taxonomia e arquitetura de domínios (quando disponível)
  ├─ context     vizinhança genômica das sequências escolhidas
  ├─ profile     uma linha comparativa por comunidade               ← bom ponto de parada
  ├─ select      opcional: aceita ou rejeita sequências usando o seu vocabulário de marcadores
  ├─ build       um HMM de perfil por grupo
  └─ scan        HMMs contra alvos + teste de separabilidade (esses grupos são famílias de verdade?)
```

**O que você recebe no final:**
- **Uma tabela de sequências**, com uma linha por sequência coletada: unidade, cluster, comunidade, taxonomia, arquitetura etc. Ela é o principal entregável.
- **Uma tabela de comunidades**, comparando todas elas: tamanho real, comprimentos, taxonomia, arquiteturas, marcadores na vizinhança e termos enriquecidos.
- **Conjuntos curados e HMMs** para os grupos que você escolher, e um **relatório de separabilidade** que mostra o quanto cada HMM reconhece o próprio grupo e não os outros.
- **Um relatório de métodos e resultados**, gerado a partir dos parâmetros realmente usados.

**Duas formas de uso:**
1. **Como ferramenta completa de linha de comando.** Os mesmos comandos de hoje rodam o fluxo inteiro.
2. **Como biblioteca Python**, por partes, dentro de um notebook ou de outra ferramenta. Dois exemplos:
   - "Tenho este FASTA de hits: me dê as comunidades."
   - "Defini estes grupos por conta própria: construa os HMMs e me diga se eles se separam."

   O segundo exemplo mostra que **testar famílias deixa de exigir uma SSN**. Qualquer agrupamento pode ser testado.

---

## 3. Regras que a ferramenta vai seguir: confirme se estão corretas do ponto de vista biológico

Estas são as regras que o redesenho trata como o núcleo do FamFlow. Hoje algumas estão implícitas ou inconsistentes no código. Depois do redesenho, ficam explícitas e testadas. **Se alguma estiver errada, precisamos saber agora.**

| # | Regra | Hoje | Proposta |
|---|---|---|---|
| R1 | **Cluster não é comunidade.** O cluster só remove redundância (sequências quase idênticas). As comunidades vêm depois, da SSN montada sobre os representantes. | Quase sempre respeitada; alguns relatórios misturam as palavras. | Separação rígida no código, nos arquivos e nos relatórios. |
| R2 | **O tamanho de uma comunidade é o número de sequências que os representantes dela representam** (a soma dos tamanhos dos clusters), e não o número de representantes. Uma comunidade com 5 representantes pode valer por 5.000 sequências. | Respeitada no `profile`. | Mantida e protegida por testes. |
| R3 | **Uma sequência que não virou representante pertence à comunidade do seu representante.** Se o representante for órfão, ela também é tratada como órfã. | **Não está definida.** Essas sequências ficam sem comunidade, o que causa erros. Por exemplo, o `select` junta todas numa comunidade falsa chamada "nan". | Registrada numa coluna própria, para não distorcer tamanhos de comunidade nem seleções (nome a confirmar, seção 6A). |
| R4 | **Todo representante está em exatamente uma comunidade, ou é órfão.** Órfãs são as candidatas mais divergentes, não lixo. | As órfãs só aparecem na tabela se você usar `--keep-orphans`. | As órfãs são sempre registradas. |
| R5 | **Se a varredura de closeness não tiver pico (poucas sequências), a unidade inteira vira uma única comunidade.** Nenhuma aresta é inventada. | Igual. | Igual, mas registrado como resultado explícito, e não só como mensagem na tela. |
| R6 | **Uma comunidade só é chamada de "família" depois que o HMM dela a separa das outras** (teste de separabilidade). Antes disso, nada é família. | Mesma ideia, mas os nomes dos HMMs e das comunidades são comparados como texto e podem não bater. | Comparação pela identidade do grupo, para o teste ser confiável. |
| R7 | **Aceitar ou rejeitar sequências só acontece na etapa de seleção, usando o seu vocabulário de marcadores.** A avaliação de vizinhança só conta marcadores e vetos; ela nunca decide. | Igual. | Igual. As regras de seleção viram uma política única, explícita e registrada (ver seção 6B). |

---

## 4. O que muda para quem usa

- **Unidades não interferem mais umas nas outras.** Hoje, se a mesma proteína aparece em duas unidades (dois domínios da mesma família, ou sementes sobrepostas), as linhas podem colidir na tabela. Agora cada sequência é identificada por *unidade + acesso*.
- **Nomes de comunidade passam a ser únicos entre unidades.** Hoje toda unidade tem o seu `Louvain_0`, `Louvain_1` etc., e as etapas que trabalham com várias unidades ao mesmo tempo (seleção, escopo de vizinhança, separabilidade) podem misturar comunidades de unidades diferentes. Agora o nome é *unidade/rótulo*.
- **O escopo "anchor" funciona mesmo quando a sua semente foi agrupada sob outro representante.** Hoje ele pode voltar vazio.
- **Dá para pedir as órfãs diretamente**, por exemplo para buscar a vizinhança genômica delas (seção 6A, "escopo orphans").
- **Resultados externos podem ser usados com segurança.** Se você rodar detecção de comunidades ou clustering em outra ferramenta, o FamFlow confere se o resultado é consistente (por exemplo, se todo representante aparece) antes de usar.
- **Os métodos de cálculo não mudam.** HMMER, mmseqs, a busca todos contra todos, a varredura de closeness e a detecção de comunidades continuam como hoje, via ROTIFER.

---

## 5. Como o código será reorganizado

### 5.1 O que muda com DDD

O redesenho segue uma abordagem chamada **DDD** (_Domain-Driven Design_, ou "projeto guiado pelo domínio"). Em termos simples, a ideia é que o código seja organizado em torno da biologia do problema, e não em torno dos programas e arquivos usados. Na prática, isso muda quatro coisas:

1. **O vocabulário do código passa a ser o vocabulário da análise.** Termos como unidade, representante, comunidade, órfã e família têm um significado único, combinado com quem entende do assunto, e o código usa exatamente esses nomes. Por isso os nomes da seção 6A importam: eles vão aparecer no código, nos arquivos e nos relatórios.
2. **As regras biológicas ficam num lugar pequeno, separado e testado.** Hoje o FamFlow é um único arquivo de cerca de 3.000 linhas, e cada comando mistura leitura de arquivos, cálculo, regras biológicas e escrita de relatório. Depois do redesenho, as regras da seção 3 (R1 a R7) ficam juntas, num núcleo pequeno, com testes próprios. Fica fácil encontrar, revisar e confirmar cada regra.
3. **As ferramentas de cálculo ficam separadas das regras.** HMMER, mmseqs, a busca todos contra todos e o alinhamento são ferramentas que fazem contas. As regras dizem o que essas contas significam. Separar as duas coisas permite trocar ou atualizar uma ferramenta sem mexer nas regras, e o contrário também.
4. **O esforço vai para onde está o diferencial.** A parte que torna o FamFlow único é a descoberta de comunidades e o teste de famílias. Essas duas partes recebem o cuidado maior. As demais, como anotação de taxonomia ou registro de execução, ficam propositalmente simples.

Com isso, as mesmas peças servem tanto para a ferramenta de linha de comando quanto para quem quiser usar o FamFlow por partes, como biblioteca Python.

### 5.2 Os módulos

O código será dividido em seis módulos principais, cada um responsável por uma parte da análise:

| Módulo | Parte da análise | Etapas de hoje | Peso no projeto |
|---|---|---|---|
| **retrieval** (busca de homólogos) | Da semente às sequências coletadas: busca por HMM, cobertura do perfil, filtros, âncora. | `seed`, `search`, `inspect`, `collect` | Apoio |
| **communities** (descoberta de comunidades) | Redução de redundância, matriz de similaridade, SSN, corte de arestas, comunidades, órfãs e o perfil de cada comunidade. | `cluster`, `matrix`, `ssn`, `profile` | **Núcleo** |
| **modeling** (modelagem de famílias) | Seleção, conjuntos curados, HMMs e teste de separabilidade. | `select`, `build`, `scan` | **Núcleo** |
| **neighborhood** (vizinhança genômica) | Busca da vizinhança e avaliação pelo vocabulário de marcadores (só contagens, nunca decisões). | `context` | Apoio |
| **annotation** (anotação) | Taxonomia e arquitetura de domínios, quando disponíveis. | `annotate` | Simples |
| **provenance** (registro de execução) | Parâmetros usados, progresso, relatórios por etapa e o relatório final de métodos e resultados. | `report`, `status` | Simples |

Os dados fluem numa só direção:

```
retrieval  →  communities  →  modeling
                   ↑               ↑
annotation ────────┘               │
neighborhood ──────────────────────┘
```

Dentro dos dois módulos do núcleo, as ferramentas de cálculo ficam numa subpasta própria, separada das regras. Além dos seis módulos, há três peças de apoio:
- **workspace**: onde ficam a tabela de sequências e os arquivos de cada etapa (a mesma estrutura de pastas de hoje);
- **api**: a porta de entrada para quem usa o FamFlow como biblioteca Python;
- **cli**: os comandos de linha de comando, que só repassam os pedidos para os módulos.

A troca para essa estrutura será gradual. A cada passo, os comandos de hoje continuam funcionando e produzindo os mesmos arquivos.

---

## 6. Decisões que dependem de você

### 6A. Nomes (vocabulário)

Queremos que o vocabulário da ferramenta seja o mesmo que você usaria para conversar sobre a análise. Os nomes ficam em inglês porque o código e os arquivos usam inglês. Para cada item, aceite, sugira outro nome ou rejeite:

| Nome proposto | O que significa | Sua decisão |
|---|---|---|
| **Membership propagation** / coluna **`represented_in`** | A regra R3: a comunidade (ou "órfã") a que uma sequência não representante pertence por meio do seu representante. | "Represented in" soa natural? Prefere outro nome? |
| **Candidate group** ("grupo candidato") | Um grupo que está sendo testado como família. Em geral é uma comunidade, mas também pode ser um grupo definido por você. | O nome está bom? (Evitamos só "group" por ser vago demais.) |
| **Selection policy** ("política de seleção") | O conjunto de regras que o `select` aplica: número mínimo de marcadores, classes de marcador exigidas, tratamento de vetos, o que fazer com sequências sem vizinhança, número mínimo de sequências por grupo. | De acordo? |
| **Community assignment** ("atribuição de comunidades") | A lista que diz em qual comunidade cada sequência está, ou se ela é órfã. Também é o que uma ferramenta externa precisaria fornecer. | De acordo? |
| **SSN outcome** ("resultado da SSN") | Tudo o que a etapa de SSN produziu para uma unidade: corte de aresta escolhido, diagnóstico da varredura, comunidades, órfãs, ou a comunidade única de último caso. | De acordo? |
| **Escopo "orphans"** | Nova opção da etapa de vizinhança: buscar a vizinhança genômica das órfãs e das sequências que elas representam. | Seria útil? |

### 6B. Classes de marcador no arquivo de vocabulário

O JSON de marcadores permite marcar cada marcador como `mandatory`, `accessory` etc. Hoje **a seleção ignora a classe**: ela só conta o número total de marcadores encontrados.

- **Opção 1 (nossa proposta):** a classe passa a valer. Exemplo: "aceitar só se houver pelo menos 1 marcador *mandatory*". Por padrão, o comportamento continua igual ao de hoje.
- **Opção 2:** tirar as classes do formato, já que hoje elas não fazem nada.

### 6C. Quais sequências entram nos conjuntos curados usados para construir os HMMs?

- **Só os representantes** (como hoje): menos redundância, menos sequências.
- **Todas as sequências coletadas:** mais dados, mais redundância.

Qual deve ser o padrão?

### 6D. Sobre quais sequências o teste de separabilidade deve rodar?

Hoje os HMMs são testados contra os mesmos representantes usados para construí-los. Isso **tende a fazer a separabilidade parecer melhor do que é**. Devemos oferecer:
- um conjunto **separado** (sequências que não entraram na construção do HMM), e/ou
- **todas as sequências coletadas**?

E qual deve ser o padrão?

### 6E. O que fazer quando a varredura de closeness fica "truncada" ou "plana"

Na escolha do corte de aresta da SSN, a curva de closeness pode ficar:
- **truncada**: ainda subindo no fim da varredura, o que indica que o melhor corte provavelmente está além da faixa testada;
- **plana**: sem pico claro, o que torna o corte escolhido quase arbitrário.

Hoje o FamFlow usa o corte automático assim mesmo e só mostra um aviso. As opções são:

1. Manter como hoje, com o resultado claramente sinalizado nos relatórios (nosso padrão).
2. **Ampliar a faixa da varredura automaticamente** quando ela vier truncada. Parece a resposta mais natural para o caso "truncada".
3. Exigir que você escolha o corte manualmente.
4. Usar uma única comunidade.

O que você gostaria que acontecesse?

---

_Ficam de fora deste documento, de propósito: detalhes de implementação, testes, coordenação de versões com o ROTIFER e o idioma dos arquivos de saída. Esses pontos estão sendo decididos à parte._
