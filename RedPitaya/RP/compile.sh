echo
swig -python $1.i
cc -O2 -fPIC -g -std=gnu99 -Wall -Werror -I/opt/redpitaya/include  -L/opt/redpitaya/lib -c  $1.c -lm -lpthread -lrp
cc -O2 -fPIC -g -std=gnu99 -Wall -Werror -I/opt/redpitaya/include  -L/opt/redpitaya/lib -c  $1_wrap.c \-I /usr/include/python2.7 -lm -lpthread -lrp
cc -O2 -fPIC -g -std=gnu99 -Wall -Werror -I/opt/redpitaya/include  -L/opt/redpitaya/lib -shared $1.o $1_wrap.o -lm -lpthread -lrp -o _$1.so
