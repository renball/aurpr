import QtQuick 2.0
import Sailfish.Silica 1.0
import neyron 1.0


Page {
    PageHeader { title: "Задание 1" }

    property double t;
    property double val;
    property double mx:0.0
    Neyron {
        id: neyron
    }



        Canvas {
            id: cvs
            anchors.fill: parent

            Timer {
            interval: 1; running: true; repeat: true
            onTriggered:{
                neyron.neyrorest()
                mx=mx+0.00005
                cvs.requestPaint()
            }
            }
            onPaint: {
                var x =0.0
                var y =0.0
                var context = getContext("2d");
                context.fillStyle = "white";
                context.fillRect(0, 0, width, height);

                context.fillStyle = 'red'
                context.beginPath()
                context.moveTo(width/20, height/3);
                while(x<mx){
                    y = neyron.neyro_dif_2(x);
                    console.log("val: " + x + " " + y)
                    context.lineTo(width /20 + x*5000, height/3-5*y);

                    x = x+0.00005
                }
                context.lineWidth = 3
                context.strokeStyle = "black"
                context.stroke();
                context.closePath();
            }
        }


//        Button {
//            onClicked: {
//                //neyron.neyro_dif()
//                t=0.0;
//                while (t<1.0){
//                    val = neyron.neyro_dif_2(t);
//                    console.log("val: " + t + " " + val)
//                    t+=0.00005;

//                }
//            }
//        }


}
