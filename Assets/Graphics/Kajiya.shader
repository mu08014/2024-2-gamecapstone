Shader "Custom/Kajiya"
{
    Properties
    {
        _DiffuseTerm("DiffuseTerm", Range(0.0, 1.0)) = 1
        _SpecularTerm("SpecularTerm", Range(0.0, 1.0)) = 1
        _SpecularPower("SpecularPower", Range(0.0, 1000.0)) = 10
        _Color ("Color", Color) = (1,1,1,1)
        _MainTex ("Texture", 2D) = "white" {}
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        //LOD 200
        
        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #include "UnityCG.cginc"
            
            float _DiffuseTerm;
            float _SpecularTerm;
            float _SpecularPower;
            fixed4 _Color;
            sampler2D _MainTex;
            
            struct appdata
            {
                float4 vertex : POSITION;
                float3 tangent : TEXCOORD1;
            };
            
            struct v2f
            {
                //float2 uv : TEXCOORD0;
                float4 pos : SV_POSITION;
                float3 tangent : TANGENT;
                float4 lightDir : TEXCOORD1;
                float4 cameraDir : TEXCOORD2;
            };
            
            v2f vert (appdata v)
            {
                v2f o;
                o.pos = UnityObjectToClipPos(v.vertex);
                o.tangent = normalize(v.tangent);
                o.lightDir = mul(unity_WorldToObject, _WorldSpaceLightPos0);
                o.cameraDir = normalize(
                    mul(unity_WorldToObject, _WorldSpaceCameraPos)
                    - v.vertex); 
                return o;
            }
            
            float4 frag (v2f i) : SV_Target
            {
                float3 diffuse = _Color * _DiffuseTerm;

                float3 halfVec = normalize(i.lightDir + i.cameraDir);
                
                float3 specular = pow( 
                    sqrt(1 - dot(i.tangent, halfVec) * dot(i.tangent, halfVec)),
                    _SpecularPower) * _SpecularTerm;
                  
                float3 color = diffuse + specular;
                
                return float4(color, 1);
            }
            ENDCG
        }
    }
    FallBack "Diffuse"
}
