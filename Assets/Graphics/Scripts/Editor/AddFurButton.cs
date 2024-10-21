using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

[CustomEditor(typeof(AddFur))]
public class AddFurButton : Editor
{
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        AddFur addFur = (AddFur)target;
        if (GUILayout.Button("Update Fur"))
        {
            addFur.initMesh();
        }
    }
}
