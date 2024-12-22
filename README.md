
## Exemplo de Workflow Usando renv:

Iniciar um projeto:

Inicialize o projeto e o renv:

renv::init()
Instalar pacotes:

Instale os pacotes que você precisa para sua análise:

install.packages("dplyr")
Realizar análise:

Carregue os pacotes e faça sua análise normalmente:

library(dplyr)
Congelar o ambiente:

Após realizar sua análise e garantir que tudo funciona corretamente, tire uma snapshot do ambiente de pacotes:

renv::snapshot()
Compartilhar o projeto:

Ao compartilhar o projeto com outros, inclua o arquivo renv.lock para que eles possam restaurar o ambiente de pacotes:

## Como configurar o ambiente do projeto caso tenha recebido este projeto

1. Instale o pacote `renv` se ainda não tiver:
   
install.packages("renv")

# Restaurar todos os pacotes

renv::restore()

# Restaurar apenas alguns pacotes
renv::restore(packages = c("dplyr", "ggplot2"))
"# TPC" 
