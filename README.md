# IPH-ECO-model

Atualizado em 20-08-2021

## Setup no Windows

1. Para de rodar o modelo IPH-ECO a partir do código fonte pela primeira vez, você precisa instalar os seguintes programas:
* Visual Studio 2019 Community, Professional ou Enterprise (https://visualstudio.microsoft.com/pt-br/downloads/)
* Ferramentas de Build do Visual Studio 2019 (https://visualstudio.microsoft.com/pt-br/downloads/)
* Intel oneAPI Base Toolkit (https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html)
* Instalar Intel oneAPI HPC Toolkit (https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit/download.html)
  
> Todos os programas listados acima são gratuitos.

2. Abrir o Visual Studio 2019

3. Abrir o projeto iph-eco-model.vfproj que está na raiz da pasta iph-eco-model

## Setup no Visual Studio 2019

1. Abrir o formulário de Propriedades do projeto que pode ser acessado em Projeto > Propriedades de iph-eco-model

2. Na opção "Debugging", inserir as seguintes informações:
* Indicar o caminho para o executável da interface gráfica (iph-eco.exe) no campo "Command" (ex. D:\IPH-ECO\iph-eco.exe)
* Indicar o caminho para que o projeto acesse as variáveis de ambiente do Windows (PATH) e o diretório com os arquivos compilados pelo projeto no campo "Environment" (ex. PATH=%PATH% "D:\iph-eco-model\x64\Debug")

3. No campo "Command Line" em "Build Events > Post-Build Events", inserir o seguinte comando:

      copy "caminho para iph-eco-model.dll" "pasta de instalação da interface gráfica"

> Ex. copy "E:\iph-eco-model\x64\Debug\iph-eco-model.dll" "E:\IPH-ECO" *

