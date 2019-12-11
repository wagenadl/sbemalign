// quickviewer.cpp

#include <stdlib.h>
#include <stdio.h>
#include <QApplication>
#include <QLabel>
#include <QScrollArea>

void usage() {
  fprintf(stderr, "Usage: quickviewer PATH [Z]\n");
  exit(1);
}

class Viewer: public QLabel {
public:
  Viewer(QString path, int z): path(path), z(z) {
    load();
  }    
  void keyPressEvent(QKeyEvent *e) {
    switch (e.key()) {
    case Qt::Key_Left:

int main(int argc, char **argv) {
  if (argc<2)
    usage();
  QString path = argv[1];
  int z = 4;
  if (argc>=3)
    z = atoi(argv[2]);

  QApplication app(argc, argv);
  
