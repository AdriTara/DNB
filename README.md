# Dynamical Network Biomarkers

This project aims to find dynamical network biomarkers from a transcriptomic expression file using a biphasic algorithm to obtain network modules in a network-interaction file.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

Python3.7 and the following libraries:
  -pandas
  -scipy
  -matplotlib
  -multiprocessing
  -datetime
that may be installed as normal python libraries, depending on the user's operating system. In debian, write on the terminal:
```
$ pip install "library"
```
Example:
```
$ pip install pandas
```
## Running the tests

This script takes as arguments:
    1-An expression file separated by tabs. Example: test-expression.csv
    2-A file determining the correspondence between the spotID and NTgenID. Example: test-probe-genid.txt
    3-A file determining the correspondence between the NtgenID and AtgenID. Example: test-orthologous.txt
    4-A network interaction file with the interactor1 and interactor2 separated by 2 tabs. Example: test-interaction-file.txt
Put the files in the same directory as the script and launch the following command in a terminal:
```
python3.7 DNBs-xx.py test-expression.csv test-probe-genid.txt test-orthologous.txt test-interaction-file.txt
```

## Authors

* **Adrián Tarazona Sánchez** - *Initial work* - [DNB](https://github.com/AdriTara/DNB)

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE.md](https://github.com/AdriTara/DNB/blob/master/LICENSE) file for details

## Acknowledgments
