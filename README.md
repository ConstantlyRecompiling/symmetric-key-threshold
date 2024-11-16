# symmetric-key-threshold


## Running MPC 
Compilation
```
/compile.py -F 256 -P 16687484526258897697932455179800985333858582303172367189654896861644253604235289 mimc [n_clients] && PLAYERS=[n_clients] Scripts/mascot.sh mimc-[n_clients] --batch-size 20
```

In separate terminal
```
python3 ExternalIO/mimc-client.py 0 2 15731825616382697382611307115081933996288172253242850321711939296329393259276
python3 ExternalIO/mimc-client.py [client id (0 indexed)] [total number of clients] [client secret value]
```